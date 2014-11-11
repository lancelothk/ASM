package edu.cwru.cbc.ASM.detect;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.detect.DataType.CpGSite;
import edu.cwru.cbc.ASM.detect.DataType.Edge;
import edu.cwru.cbc.ASM.detect.DataType.MappedRead;
import edu.cwru.cbc.ASM.detect.DataType.Vertex;
import edu.cwru.cbc.ASM.detect.Utils.MappedReadFileLineProcessor;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by ke on 2/19/14.
 * ASM detection
 * Note:    1. mapped reads start pos is 1-based, end pos is 0-based.
 */
public class DetectASM implements Callable<String> {
    private static final int MIN_READ_CPG = 2;
    private File intervalFile;
    private long initPos;
    private BufferedWriter groupWriter;
    private BufferedWriter group2Writer;

    public DetectASM(File intervalFile, long initPos, BufferedWriter writer, BufferedWriter group2Writer) {
        this.intervalFile = intervalFile;
        this.initPos = initPos;
        this.groupWriter = writer;
        this.group2Writer = group2Writer;
    }

    public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
        long start = System.currentTimeMillis();
        // call on single file
//        detectASM.execute("/home/kehu/ASM_result/chr20-56897421-56898208.reads", 56897421);
//        detectASM.execute("", 1);
//		detectASM.execute(new File("/home/lancelothk/chr20_test/chr20-56895353-56895567"), 56895353);
        // test reads setting
        String pathName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/chrTest2-1-6";
        String summaryFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.summary";
        String groupResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.groupResult";
        String group2ResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.group2Result";

        // TODO mkdir if not exist
        String cellLine = "i90";
        String replicate = "r1";
        String name = "_atLeastTwo_large";

//        String pathName = String.format("/home/kehu/experiments/ASM/result_%s_%s/intervals%s", cellLine, replicate,
//                                        name);
//
//        String summaryFileName = String.format(
//                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_summary%3$s", cellLine, replicate,
//                name);
//        String groupResultFileName = String.format(
//                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_groups%3$s", cellLine, replicate,
//                name);
//        String group2ResultFileName = String.format(
//                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_group2%3$s", cellLine, replicate,
//                name);

        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
        summaryWriter.write("name\tlength\treadCount\tCpGCount\tGroupCount\n");
        BufferedWriter groupWriter = new BufferedWriter(new FileWriter(groupResultFileName));
        BufferedWriter group2Writer = new BufferedWriter(new FileWriter(group2ResultFileName));
        group2Writer.write("name\tlength\treadCount\tCpGCount\tGroupCount\t1stGroupSize\t2ndGroupSize\n");
        File path = new File(pathName);

        ExecutorService executor = Executors.newFixedThreadPool(4);
        List<Future<String>> futureList = new ArrayList<>();
        if (path.isDirectory()) {
            for (File file : path.listFiles()) {
                if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
                        !file.getName().endsWith("intervalSummary")) {
                    String[] items = file.getName().split("-");
                    Future<String> future = executor.submit(
                            new DetectASM(file, Long.parseLong(items[1]), groupWriter, group2Writer));
                    futureList.add(future);
                }
            }
        } else {
            String[] items = path.getName().split("-");
            Future<String> future = executor.submit(
                    new DetectASM(path, Long.parseLong(items[1]), groupWriter, group2Writer));
            futureList.add(future);
        }

        for (Future<String> stringFuture : futureList) {
            summaryWriter.write(stringFuture.get());
            //System.out.println(stringFuture.get());
        }
        executor.shutdown();

        summaryWriter.close();
        groupWriter.close();
        group2Writer.close();
        System.out.println(System.currentTimeMillis() - start + "ms");
    }


    @Override
    public String call() throws Exception {
        StringBuilder asm_result = new StringBuilder(intervalFile.getName() + "\n");
        StringBuilder summary = new StringBuilder(intervalFile.getName());
        String reference = readRef(intervalFile);
        List<CpGSite> cpgList = extractCpGSite(reference, initPos);
        List<MappedRead> mappedReadList = Files.asCharSource(intervalFile, Charsets.UTF_8).readLines(
                new MappedReadFileLineProcessor());
        associateReadWithCpG(cpgList, mappedReadList);
        String[] items = intervalFile.getName().split("-");
        if (items.length != 3) {
            throw new RuntimeException("illegal interval file name format!");
        }
        summary.append("\t" + (Integer.parseInt(items[2]) - Integer.parseInt(items[1]) + 1));
        summary.append("\t" + mappedReadList.size());
        summary.append("\t" + cpgList.size());

        List<Edge> edgeList = new ArrayList<>();
        Map<Long, Vertex> vertexMap = new HashMap<>();

        constructGraph(vertexMap, edgeList, mappedReadList);
        System.out.println(
                intervalFile.getName() + "\t" + "vertex number:" + vertexMap.keySet().size() + "\t" + "edge number:" +
                        edgeList.size());
        getClusters(asm_result, edgeList, vertexMap);

        int groupCount = 0;

        for (Vertex vertex : vertexMap.values()) {
            asm_result.append("group: " + groupCount++ + "\n");
            Collections.sort(vertex.getMappedReadList(),
                             (MappedRead read1, MappedRead read2) -> (int) (read1.getStart() - read2.getStart()));
            for (MappedRead mappedRead : vertex.getMappedReadList()) {
                asm_result.append(mappedRead.getId() + ",");
            }
            asm_result.append("\n");
            asm_result.append("Covered CpG sites:\n");
            List<CpGSite> cpGSiteList = new ArrayList<>(vertex.getCoveredCpGSites());
            Collections.sort(cpGSiteList, (CpGSite c1, CpGSite c2) -> (int) (c1.getPos() - c2.getPos()));
            for (CpGSite cpGSite : cpGSiteList) {
                double methylCount = 0, totalCount = 0;
                for (MappedRead mappedRead : vertex.getMappedReadList()) {
                    if (mappedRead.getFirstCpG().getPos() <= cpGSite.getPos() &&
                            mappedRead.getLastCpG().getPos() >= cpGSite.getPos()) {
                        if (mappedRead.getCpGMethylStatus(cpGSite.getPos())) {
                            methylCount++;
                        }
                        totalCount++;
                    }
                }
                asm_result.append(String.format("%d-%.2f,", cpGSite.getPos(), methylCount / totalCount));
            }
            asm_result.append("\n");
        }
        asm_result.append("Number of groups:\t" + groupCount + "\n");
        summary.append("\t" + groupCount + "\n");
        String intervalSummary = summary.toString();
        if (vertexMap.values().size() == 2) {
            group2Writer.write(intervalSummary + "\t");
            for (Vertex vertex : vertexMap.values()) {
                group2Writer.write(vertex.getMappedReadList().size() + "\t");
            }
            group2Writer.write("\n");
        }
        groupWriter.write(asm_result.toString());

        Map<Long, Integer> cpgPosMap = cpgList.stream().collect(Collectors.toMap(CpGSite::getPos, cpgSite -> 0));
        for (Vertex vertex : vertexMap.values()) {
            for (CpGSite cpGSite : vertex.getCoveredCpGSites()) {
                if (cpgPosMap.containsKey(cpGSite.getPos())) {
                    cpgPosMap.put(cpGSite.getPos(), cpgPosMap.get(cpGSite.getPos()) + 1);
                }
            }
        }

        SortedSet<Long> sortedPosSet = new TreeSet<>(cpgPosMap.keySet());
        for (Long aLong : sortedPosSet) {
            groupWriter.write(String.format("pos: %d\tcount: %d\n", aLong, cpgPosMap.get(aLong)));
        }
        return summary.toString();
    }

    private void getClusters(StringBuilder asm_result, List<Edge> edgeList,
                             Map<Long, Vertex> vertexMap) throws IOException {
        long start = System.currentTimeMillis();
        // count how many tie situation occurs. For analysis use
        int tieWeightCounter = 0, tieIdCountCounter = 0, tieMethylPolarity = 0;
        // merge vertexes connected by positive weight edge
        while (true) {
            // if edgeList is empty
            if (edgeList.size() == 0) {
                break;
            }
            List<Edge> maxEdgeList = getMaxFromList(edgeList, Edge::getWeight);
            // if max weight <= 0, stop merge.
            if (maxEdgeList.get(0).getWeight() <= 0) {
                break;
            } else {
                // since edgeList is not empty, getMaxFromList always returns at least one element
                if (maxEdgeList.size() != 1) {
                    // multiple equal max weight edges.
                    tieWeightCounter++;
                    // use methyl polarity to group similar vertex first
                    // sort edge by abs of methyl polarity. Larger first
//                    maxEdgeList = getMaxFromList(maxEdgeList, Edge::getMethylPolarityAbs);
//                    if (maxEdgeList.size() != 1) {
//                        tieMethylPolarity++;

                    // First pick edge not connect to two clusters.
                    // sort edge by id count. Smaller first
                    maxEdgeList = getMinFromList(maxEdgeList, Edge::getIdCount);
                    if (maxEdgeList.size() != 1) {
                        tieIdCountCounter++;
                    }
//                    }
                }
                mergeVertex(maxEdgeList.get(0), edgeList, vertexMap);
            }
        }
        asm_result.append(String.format("tied weight counter:%d\n", tieWeightCounter));
        asm_result.append(String.format("tied methyl polarity:%d\n", tieMethylPolarity));
        asm_result.append(String.format("tied id counter:%d\n", tieIdCountCounter));
    }

    private void mergeVertex(Edge edge, List<Edge> edgeList, Map<Long, Vertex> vertexMap) {
        Vertex left = edge.getLeft();
        Vertex right = edge.getRight();
        // merge right to left vertex
        left.addMappedRead(right.getMappedReadList());
        left.addInnerWeight(edge.getWeight());
        left.addCpGSites(right.getCoveredCpGSites());
        //remove this edge and right vertex from graph
        edge.removeFromVertex();
        edgeList.remove(edge);
        vertexMap.remove(right.getId());
        // update right vertex adj edges
        for (Edge edgeR : right.getAdjEdges()) {
            // update new vertex
            if (!edgeR.replaceVertex(right, left)) {
                throw new RuntimeException("fail to update new vertex in adj edge!");
            }
        }
        // update left vertex with new edges
        for (Edge adjEdge : right.getAdjEdges()) {
            left.addEdge(adjEdge);
        }
        // update edges of vertex which connect both left and right
        updateAndRemoveDupEdge(edgeList);
    }

    private void updateAndRemoveDupEdge(List<Edge> edgeList) {
        Map<String, Edge> edgeMap = new HashMap<>();
        Iterator<Edge> edgeIterator = edgeList.iterator();
        while (edgeIterator.hasNext()) {
            Edge edgeB = edgeIterator.next();
            // duplicate edge
            if (edgeMap.containsKey(edgeB.getUniqueId())) {
                Edge edgeA = edgeMap.get(edgeB.getUniqueId());
                // keep edgeA, remove edgeB.
                edgeA.setWeight(edgeA.getWeight() + edgeB.getWeight());
                edgeB.removeFromVertex();
                edgeIterator.remove(); // remove from edgelist
            } else {
                edgeMap.put(edgeB.getUniqueId(), edgeB);
            }
        }
    }

    /**
     * select items with max property value in the list.
     *
     * @param itemList
     * @param getProperty
     * @param <T>
     * @return list contains items with max weight
     */
    private <T> List<T> getMaxFromList(List<T> itemList, Function<T, Integer> getProperty) {
        List<T> maxItemList = new ArrayList<>();
        int maxValue = Integer.MIN_VALUE;
        for (T item : itemList) {
            if (getProperty.apply(item) > maxValue) {
                maxItemList.clear();
                maxItemList.add(item);
                maxValue = getProperty.apply(item);
            } else if (getProperty.apply(item) == maxValue) {
                maxItemList.add(item);
            }
        }
        return maxItemList;
    }

    /**
     * select items with min property value in the list.
     *
     * @param itemList
     * @param getProperty
     * @param <T>
     * @return list contains items with min weight
     */
    private <T> List<T> getMinFromList(List<T> itemList, Function<T, Integer> getProperty) {
        List<T> minItemList = new ArrayList<>();
        int minValue = Integer.MAX_VALUE;
        for (T item : itemList) {
            if (getProperty.apply(item) < minValue) {
                minItemList.clear();
                minItemList.add(item);
                minValue = getProperty.apply(item);
            } else if (getProperty.apply(item) == minValue) {
                minItemList.add(item);
            }
        }
        return minItemList;
    }

    private void constructGraph(Map<Long, Vertex> vertexMap, List<Edge> edgeList, List<MappedRead> mappedReadList) {
        // visit all possible edges once.
        for (int i = 0; i < mappedReadList.size(); i++) {
            // only use reads which contains at least two Cpg sites
            if (mappedReadList.get(i).getCpgList().size() >= MIN_READ_CPG) {
                for (int j = i + 1; j < mappedReadList.size(); j++) {
                    int score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
                    if (score != Integer.MIN_VALUE) {
                        long idI = mappedReadList.get(i).getId(), idJ = mappedReadList.get(j).getId();
                        if (!vertexMap.containsKey(idI)) {
                            vertexMap.put(idI, new Vertex(mappedReadList.get(i)));
                        }
                        if (!vertexMap.containsKey(idJ)) {
                            vertexMap.put(idJ, new Vertex(mappedReadList.get(j)));
                        }
                        edgeList.add(new Edge(vertexMap.get(idI), vertexMap.get(idJ), score));
                    }
                }
            }
        }
    }

    private int checkCompatible(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            int score = 0;
            for (CpGSite cpgA : readA.getCpgList()) {
                for (CpGSite cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos()) {
                        if (cpgA.isMethylated() != cpgB.isMethylated()) {
                            score--;
                        } else {
                            score++;
                        }
                    }
                }
            }
            return score;
        } else {
            // don't have overlapped CpG, non-connected
            return Integer.MIN_VALUE;
        }
    }

    private String readRef(File intervalFile) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(intervalFile));
        String line = bufferedReader.readLine();
        String[] items = line.split("\t");
        if (items.length != 2) {
            throw new RuntimeException("invalid reference line in interval read file!");
        } else {
            return items[1].toUpperCase();
        }
    }

    private void associateReadWithCpG(List<CpGSite> cpgList, List<MappedRead> mappedReadList) {
        for (MappedRead read : mappedReadList) {
            for (CpGSite cpg : cpgList) {
                if (cpg.getPos() >= read.getStart() && cpg.getPos() + 1 <= read.getEnd()) {
                    cpg.addAssociatedRead(read);
                    read.addCpG(new CpGSite(cpg.getPos(), read.getCpGMethylStatus(cpg.getPos())));
                }
            }
        }
    }

    private List<CpGSite> extractCpGSite(String reference, long initPos) {
        reference = reference.replace(" ", "");
        List<CpGSite> cpgList = new ArrayList<>();
        for (int i = 0; i < reference.length() - 1; i++) {
            if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
                cpgList.add(new CpGSite(initPos + i));
            }
        }
        return cpgList;
    }

}