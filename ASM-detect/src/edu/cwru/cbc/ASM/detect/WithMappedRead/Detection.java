package edu.cwru.cbc.ASM.detect.WithMappedRead;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.ASMGraph;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.ClusterRefCpG;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.GroupResult;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.Vertex;
import edu.cwru.cbc.ASM.visualization.ReadsVisualization;
import org.apache.commons.cli.*;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;

/**
 * Created by kehu on 11/12/14.
 * ASM Detection with whole read info.
 */
public class Detection implements Callable<String> {
    private File inputFile;
    private String chr;
    private int startPos;
    private int endPos;
    private int min_interval_cpg;
    private double fisher_p_threshold;
    private double region_threshold;
    private boolean useRegionP;

    /**
     * Detection constructor.
     *
     * @param inputFile          A file contains reads in input region or a folder contains all input region files. File name should be in format:chr-start-end
     * @param min_interval_cpg   Minimum number of CpGs in the interval. If under this threshold(too small), won't compute the error probability.
     * @param fisher_p_threshold Significant fisher P value threshold.
     * @param region_threshold   Significant region P value threshold.
     * @param useRegionP         If use the region P value as result filtering criteria or use region percent.
     */
    public Detection(File inputFile, int min_interval_cpg, double fisher_p_threshold,
                     double region_threshold, boolean useRegionP) {
        this.fisher_p_threshold = fisher_p_threshold;
        this.region_threshold = region_threshold;
        this.inputFile = inputFile;
        this.min_interval_cpg = min_interval_cpg;
        this.useRegionP = useRegionP;
    }

    public static void main(
            String[] args) throws ParseException, IOException, ExecutionException, InterruptedException {
        long start = System.currentTimeMillis();

        Options options = new Options();
        options.addOption("i", true, "Input intervals folder or interval file name");
        options.addOption("s", true, "Summary File");
        options.addOption("mic", true, "Minimum interval cpg number");
        options.addOption("mrc", true, "Minimum read CpG number");
        options.addOption("f", true, "Fisher exact test P threshold");
        options.addOption("m", true, "Mode: Percent(per) or RegionP(rp)");
        options.addOption("r", true, "Region threshold");
        options.addOption("t", false, "Thread number to execute the program.");

        CommandLineParser parser = new BasicParser();
        CommandLine cmd = parser.parse(options, args);

        String inputPathName = cmd.getOptionValue("i");
        String summaryFileName = cmd.getOptionValue("s");
        int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
        int min_read_cpg = Integer.valueOf(cmd.getOptionValue("mrc"));
        double fisher_p_threshold = Double.valueOf(cmd.getOptionValue("f"));
        double region_threshold;
        boolean useRegionP;
        if (cmd.getOptionValue("m").equals("per")) {
            region_threshold = Double.valueOf(cmd.getOptionValue("r"));
            useRegionP = false;
        } else if (cmd.getOptionValue("m").equals("rp")) {
            region_threshold = Double.valueOf(cmd.getOptionValue("r"));
            useRegionP = true;
        } else {
            throw new RuntimeException("invalid mode!");
        }
        int threadNumber = Integer.valueOf(cmd.getOptionValue("t", "6"));
        List<String> resultList = execute(inputPathName, threadNumber, min_interval_cpg, min_read_cpg,
                                          fisher_p_threshold, region_threshold, useRegionP);
        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
        summaryWriter.write(
                "chr\tstartPos\tendPos\tlength\tvertex number\tedge number\treadCount\t#CpGSite\t#CpGSiteInClusters\tGroupCount\t" +
                        "avgGroupPerCpG\tMECSum\tCpGSum\tNormMEC\terrorProbability\t#CpGwithFisherP<=" +
                        fisher_p_threshold + "\t" + (useRegionP ? "regionP" : "RegionPercent") + "=" +
                        region_threshold + "\tgroup1\tgroup2\tlabel\n");

        for (String result : resultList) {
            summaryWriter.write(result);
        }
        summaryWriter.close();
        System.out.println(System.currentTimeMillis() - start + "ms");
    }

    private static List<String> execute(String inputName, int threadNumber, int min_interval_cpg, int min_read_cpg,
                                        double fisher_p_threshold, double region_threshold,
                                        boolean useRegionP) throws ExecutionException, InterruptedException {
        File inputFile = new File(inputName);
        List<String> resultList = new ArrayList<>();
        ExecutorService executor = Executors.newFixedThreadPool(threadNumber);
        List<Future<String>> futureList = new ArrayList<>();
        if (inputFile.isDirectory()) {
            File[] files = inputFile.listFiles();
            if (files == null) {
                throw new RuntimeException("Empty folder!");
            } else {
                for (File file : files) {
                    try {
                        if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
                                !file.getName().endsWith("test") && !file.getName().endsWith("~") &&
                                !file.getName().endsWith("detected")) {
                            Future<String> future = executor.submit(
                                    new Detection(file, min_interval_cpg, fisher_p_threshold,
                                                  region_threshold, useRegionP));
                            futureList.add(future);
                        }
                    } catch (Exception e) {
                        throw new RuntimeException("Problem File name: " + file.getAbsolutePath() + "\n", e);
                    }
                }
            }
        } else {
            try {
                Future<String> future = executor.submit(
                        new Detection(inputFile, min_interval_cpg, fisher_p_threshold, region_threshold,
                                      useRegionP));
                futureList.add(future);
            } catch (Exception e) {
                throw new RuntimeException("Problem File name:" + inputFile.getAbsolutePath(), e);
            }
        }
        for (Future<String> stringFuture : futureList) {
            resultList.add(stringFuture.get());
        }
        executor.shutdown();
        return resultList;
    }

    private void extractIntervalPosition(File inputFile) {
        String[] items = inputFile.getName().split("-");
        if (items.length != 3) {
            throw new RuntimeException("invalid input file name format!\t" + inputFile.getName());
        }
        chr = items[0];
        startPos = Integer.parseInt(items[1]);
        endPos = Integer.parseInt(items[2]);
    }

    public String call() throws Exception {
        extractIntervalPosition(inputFile);

        // load input
        String reference = DetectionUtils.readRefFromIntervalFile(inputFile);
        List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
        // filter out reads which only cover 1 or no CpG sites
        List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8).readLines(
                new MappedReadLineProcessor(refCpGList));

        // align read
        ReadsVisualization.alignReads(mappedReadList, reference, inputFile.getAbsolutePath() + ".aligned");

        // construct graph
        ASMGraph graph = new ASMGraph(mappedReadList);

        // clustering
        graph.cluster();

        List<RefCpG> twoClusterRefCpGList = getTwoClustersRefCpG(refCpGList, graph.getClusterRefCpGMap());
        // get fisher test P values for each refCpG in clusters.
        double regionP;
        if (fisherTest(graph, twoClusterRefCpGList)) {
            // use product of 1-pi
//            regionP = 1;
//            for (RefCpG refCpG : twoClusterRefCpGList) {
//                regionP *= (1 - refCpG.getP_value());
//            }
// 			  regionP = 1 - regionP;

            // use min P -- sidak's combination test
//			double minP = Double.MAX_VALUE;
//			for (RefCpG refCpG : twoClusterRefCpGList) {
//				if (minP > refCpG.getP_value()){
//					minP = refCpG.getP_value();
//				}
//			}
//			regionP = 1- Math.pow(1-minP, twoClusterRefCpGList.size());

            // use fisher's combination test
            regionP = 0;
            for (RefCpG refCpG : twoClusterRefCpGList) {
                regionP += Math.log(refCpG.getP_value());
            }
            ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(2 * twoClusterRefCpGList.size());
            regionP = 1 - chiSquaredDistribution.cumulativeProbability(-2 * regionP);

            // use Stouffer's combination test
//            double z=0;
//            NormalDistribution stdNorm = new NormalDistribution(0,1);
//            for (RefCpG refCpG : twoClusterRefCpGList) {
//                z += stdNorm.inverseCumulativeProbability(1 - refCpG.getP_value());
//            }
//            z /=Math.sqrt(twoClusterRefCpGList.size());
//            regionP = 1 - stdNorm.cumulativeProbability(z);
        } else {
            // give -1 if only one cluster
            regionP = -1;
        }

        double avgGroupCpGCoverage = writeGroupResult(refCpGList, graph);

        return buildSummary(endPos - startPos + 1, refCpGList, mappedReadList, graph, avgGroupCpGCoverage,
                            twoClusterRefCpGList, regionP);
    }

    private double writeGroupResult(List<RefCpG> refCpGList, ASMGraph graph) throws IOException {
        BufferedWriter groupResultWriter = new BufferedWriter(
                new FileWriter(inputFile.getAbsolutePath() + ".detected"));

        groupResultWriter.write(inputFile.getName() + "\n");
        groupResultWriter.write(String.format("tied weight counter:%d\n", graph.getTieWeightCounter()));
        groupResultWriter.write(String.format("tied id counter:%d\n", graph.getTieIdCountCounter()));

        double sum = printRefCpGGroupCoverage(refCpGList, graph, groupResultWriter);

        List<GroupResult> groupResultList = writeAlignedResult(refCpGList, graph, groupResultWriter);
        writePvalues(refCpGList, groupResultWriter);

        writeGroupReadList(groupResultWriter, groupResultList);

        groupResultWriter.close();

        return sum / graph.getClusterRefCpGMap().size();
    }

    private double printRefCpGGroupCoverage(List<RefCpG> refCpGList, ASMGraph graph,
                                            BufferedWriter groupResultWriter) throws IOException {
        // print refCpG's group coverage
        double sum = 0;
        for (RefCpG refCpG : refCpGList) {
            if (graph.getClusterRefCpGMap().containsKey(refCpG.getPos())) {
                groupResultWriter.write(String.format("pos: %d\tcount: %d\n", refCpG.getPos(),
                                                      graph.getClusterRefCpGMap().get(
                                                              refCpG.getPos()).getClusterCount()));
                sum += graph.getClusterRefCpGMap().get(refCpG.getPos()).getClusterCount();
            }
        }
        return sum;
    }

    private List<GroupResult> writeAlignedResult(List<RefCpG> refCpGList, ASMGraph graph,
                                                 BufferedWriter groupResultWriter) throws IOException {
        List<GroupResult> groupResultList = graph.getClusterResult().values().stream().map(
                vertex -> new GroupResult(new ArrayList<>(vertex.getRefCpGMap().values()), vertex.getMappedReadList(),
                                          vertex.getMECScore())).collect(Collectors.toList());

        // start to write aligned result
        groupResultWriter.write("Aligned result:\n");
        groupResultList.sort(GroupResult::compareTo);

        for (int i = 0; i < groupResultList.size(); i++) {
            groupResultWriter.write(i + ":\t");
            for (RefCpG refCpG : refCpGList) {
                Map<Integer, RefCpG> refCpGMap = groupResultList.get(i).getRefCpGList().stream().collect(
                        Collectors.toMap(RefCpG::getPos, r -> r));
                if (refCpGMap.containsKey(refCpG.getPos())) {
                    groupResultWriter.write(Strings.padEnd(
                            String.format("%.2f(%d)", refCpGMap.get(refCpG.getPos()).getMethylLevel(),
                                          refCpGMap.get(refCpG.getPos()).getCoveredCount()), 10, ' '));
                } else {
                    // fill the gap
                    groupResultWriter.write(Strings.repeat(" ", 10));
                }
            }
            groupResultWriter.write("\n");
        }
        return groupResultList;
    }

    private void writePvalues(List<RefCpG> refCpGList, BufferedWriter groupResultWriter) throws IOException {
        groupResultWriter.write("P:\t");
        for (RefCpG refCpG : refCpGList) {
            if (refCpG.getP_value() != -1) {
                groupResultWriter.write(Strings.padEnd(String.format("%.3f", refCpG.getP_value()), 10, ' '));
            } else {
                groupResultWriter.write(Strings.repeat(" ", 10));

            }
        }
        groupResultWriter.write("\n");
    }

    private void writeGroupReadList(BufferedWriter groupResultWriter,
                                    List<GroupResult> groupResultList) throws IOException {
        groupResultWriter.write("Group's read list:\n");
        for (int i = 0; i < groupResultList.size(); i++) {
            groupResultWriter.write(
                    "Group:" + i + "\t" + "size:" + groupResultList.get(i).getMappedReadList().size() + "\tMEC:" +
                            groupResultList.get(i).getMec() + "\n");
            for (MappedRead mappedRead : groupResultList.get(i).getMappedReadList()) {
                groupResultWriter.write(mappedRead.getId() + ",");
            }
            groupResultWriter.write("\n");
        }
    }

    private int getFisherTestCount(List<RefCpG> twoClusterRefCpGList) {
        int testCount = 0;
        for (RefCpG refCpG : twoClusterRefCpGList) {
            if (refCpG.getP_value() <= this.fisher_p_threshold) {
                testCount++;
            }
        }
        return testCount;
    }

    private boolean fisherTest(ASMGraph graph, List<RefCpG> twoClusterRefCpGList) {
        if (graph.getClusterResult().size() != 1) { // cluster size == 2 or more
            for (RefCpG refCpG : twoClusterRefCpGList) {
                assert refCpG.getCpGCoverage() == 2;
                int j = 0;
                int[][] matrix = new int[2][2];
                for (Vertex vertex : graph.getClusterResult().values()) {
                    if (vertex.getRefCpGMap().containsKey(refCpG.getPos())) {
                        matrix[j][0] = vertex.getRefCpGMap().get(refCpG.getPos()).getMethylCount();
                        matrix[j][1] = vertex.getRefCpGMap().get(refCpG.getPos()).getNonMethylCount();
                        j++;
                    }
                }
                double fisher_P = FisherExactTest.fishersExactTest(matrix[0][0], matrix[0][1], matrix[1][0],
                                                                   matrix[1][1])[0];  // [0] is two tail test.
                refCpG.setP_value(fisher_P);
            }
            return true;
        } else {
            return false;
        }
    }

    private String buildSummary(int length, List<RefCpG> refCpGList, List<MappedRead> mappedReadList, ASMGraph graph,
                                double avgGroupCpGCoverage, List<RefCpG> twoClusterRefCpGList, double regionP) {
        int testCount = getFisherTestCount(twoClusterRefCpGList);
        StringBuilder sb = new StringBuilder(
                String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%f\t%f\t%d\t%.5f\t", chr, startPos,
                              endPos, length, graph.getOriginalVertexCount(), graph.getOriginalEdgeCount(),
                              mappedReadList.size(), refCpGList.size(), twoClusterRefCpGList.size(),
                              graph.getClusterResult().size(), avgGroupCpGCoverage, graph.getMECSum(),
                              graph.getCpGSum(), graph.getNormMECSum(),
                              calcErrorProbability(graph.getClusterResult().values(), twoClusterRefCpGList), testCount,
                              useRegionP ? regionP : testCount / (double) refCpGList.size()));
        switch (graph.getClusterResult().values().size()) {
            case 1:
                for (Vertex vertex : graph.getClusterResult().values()) {
                    sb.append(vertex.getMappedReadList().size()).append("\t");
                }
                sb.append("NULL\t");
                break;
            case 2:
                for (Vertex vertex : graph.getClusterResult().values()) {
                    sb.append(vertex.getMappedReadList().size()).append("\t");
                }
                break;
            default:
                throw new RuntimeException("more than 2 clusters in result!");

        }
        if (useRegionP) {
            if (regionP <= this.region_threshold && regionP >= 0) {
                sb.append('+');
            } else {
                sb.append('-');
            }
        } else {
            if (testCount / (double) refCpGList.size() >= this.region_threshold) {
                sb.append('+');
            } else {
                sb.append('-');
            }
        }

        sb.append("\n");
        return sb.toString();
    }

    private double calcErrorProbability(Collection<Vertex> clusterResult, List<RefCpG> twoClusterRefCpGList) {
        double p;
        switch (clusterResult.size()) {
            case 1:
                p = -1; // one cluster case
                break;
            case 2:
                List<Vertex> clusterResultList = clusterResult.stream().collect(Collectors.toList());
                clusterResultList.sort(
                        (Vertex v1, Vertex v2) -> v1.getMappedReadList().size() - v2.getMappedReadList().size());
                Vertex minorityCluster = clusterResultList.get(0);
                Vertex majorityCluster = clusterResultList.get(1);

                // overlapped CpG < GlobalParameter.MIN_INTERVAL_CPG
                if (twoClusterRefCpGList.size() < this.min_interval_cpg) {
                    p = -3;
                } else {
                    p = EPCaluculation_average_combine(minorityCluster, majorityCluster, twoClusterRefCpGList);
                }
                break;
            default:
                p = -2; // more than 2 clusters
        }
        return p;
    }

    private List<RefCpG> getTwoClustersRefCpG(List<RefCpG> refCpGList, Map<Integer, ClusterRefCpG> clusterRefCpGMap) {
        return refCpGList.stream().filter(refCpG -> clusterRefCpGMap.containsKey(refCpG.getPos())).filter(
                refCpG -> clusterRefCpGMap.get(refCpG.getPos()).getClusterCount() == 2).collect(Collectors.toList());
    }

    private double EPCaluculation_average_combine(Vertex minorityCluster, Vertex majorityCluster,
                                                  List<RefCpG> twoClusterRefCpGList) {
        double p = 0;
        for (RefCpG refCpG : twoClusterRefCpGList) {
            double majorMethylLevel = majorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylLevel();
            long comb = CombinatoricsUtils.binomialCoefficient(
                    minorityCluster.getRefCpGMap().get(refCpG.getPos()).getCoveredCount(),
                    minorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylCount());
            p += comb *
                    FastMath.pow(majorMethylLevel,
                                 minorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylCount()) *
                    FastMath.pow(1 - majorMethylLevel,
                                 minorityCluster.getRefCpGMap().get(refCpG.getPos()).getNonMethylCount());
        }
        p /= twoClusterRefCpGList.size();
        return p;
    }
}
