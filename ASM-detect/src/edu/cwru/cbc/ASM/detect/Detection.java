package edu.cwru.cbc.ASM.detect;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.IO.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import edu.cwru.cbc.ASM.detect.DataType.*;
import edu.cwru.cbc.ASM.visualization.ReadsVisualization;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
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
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils.extractCpGSite;

/**
 * Created by lancelothk on 6/8/15.
 * ASM Detection with whole read info.
 */
public class Detection implements Callable<IntervalDetectionSummary> {
	private final int ALIGN_COL_SIZE = 12;
	private File inputFile;
	private String chr;
	private int startPos;
	private int endPos;
	private int min_interval_cpg;

	/**
	 * Detection constructor.
	 *
	 * @param inputFile        A file contains reads in input region or a folder contains all input region files. File name should be in format:chr-start-end
	 * @param min_interval_cpg Minimum number of CpGs in the interval. If under this threshold(too small), won't compute the error probability.
	 */
	public Detection(File inputFile, int min_interval_cpg) {
		this.inputFile = inputFile;
		this.min_interval_cpg = min_interval_cpg;
	}

	public IntervalDetectionSummary call() throws Exception {
		extractIntervalPosition(inputFile);

		// load input
		String reference = DetectionUtils.readRefFromIntervalFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		// filter out reads which only cover 1 or no CpG sites
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor(refCpGList));

		// construct graph
		ASMGraph graph = new ASMGraph(mappedReadList);

		// clustering
		graph.cluster();

		List<RefCpG> twoClusterRefCpGList = getTwoClustersRefCpG(refCpGList, graph.getClusterRefCpGMap());

		double regionP;
		if (twoClusterRefCpGList.size() < min_interval_cpg) {
			// give 3 if interval contain less #cpg than min_interval_cpg
			regionP = 3;
		} else if (fisherTest(graph, twoClusterRefCpGList)) {
			// get fisher test P values for each refCpG in clusters.
			regionP = calcRegionP_StoufferComb(twoClusterRefCpGList);

			// update start/end position for detected AMR region. Excluding single cluster CpG in the boundary.
			twoClusterRefCpGList.sort((r1, r2) -> r1.getPos() - r2.getPos());
			startPos = twoClusterRefCpGList.get(0).getPos();
			endPos = twoClusterRefCpGList.get(twoClusterRefCpGList.size() - 1).getPos() + 1;
		} else {
			// give 2 if only one cluster
			regionP = 2;
		}

		List<GroupResult> groupResultList = writeGroupResult(refCpGList, graph);
		if (graph.getClusterResult().values().size() != 2) {
			throw new RuntimeException("clusters number is not 2!");
		}
		List<List<MappedRead>> readGroups = graph.getClusterResult().values().stream().map(
				Vertex::getMappedReadList).collect(
				Collectors.toList());
		ReadsVisualization.alignReadsIntoGroups(readGroups, reference, inputFile.getAbsolutePath() + ".groups.aligned");
		return new IntervalDetectionSummary(regionP, chr.replace("chr", ""), startPos, endPos, endPos - startPos + 1,
				graph.getOriginalEdgeCount(), mappedReadList.size(), refCpGList.size(),
				twoClusterRefCpGList.size(), graph.getClusterResult().size(), graph.getCpGSum(),
				graph.getMECSum(), graph.getNormMECSum(),
				calcErrorProbability(graph.getClusterResult().values(), twoClusterRefCpGList), regionP,
				Iterables.get(graph.getClusterResult().values(), 0).getMappedReadList().size(),
				Iterables.get(graph.getClusterResult().values(), 1).getMappedReadList().size(),
				groupResultList.get(0).getAvgMethylLevel(),
				groupResultList.get(1).getAvgMethylLevel(),
				"<label>");
	}

	private void extractIntervalPosition(File inputFile) {
		String[] items = inputFile.getName().replace(Constant.MAPPEDREADS_EXTENSION, "").split("-");
		if (items.length != 3) {
			throw new RuntimeException("invalid input file name format!\t" + inputFile.getName());
		}
		chr = items[0];
		startPos = Integer.parseInt(items[1]);
		endPos = Integer.parseInt(items[2]);
	}

	private List<RefCpG> getTwoClustersRefCpG(List<RefCpG> refCpGList, Map<Integer, ClusterRefCpG> clusterRefCpGMap) {
		return refCpGList.stream()
				.filter(refCpG -> clusterRefCpGMap.containsKey(refCpG.getPos()))
				.filter(refCpG -> clusterRefCpGMap.get(refCpG.getPos()).getClusterCount() == 2)
				.collect(Collectors.toList());
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
				if (fisher_P >= 0 && fisher_P <= 1) {
					refCpG.setP_value(fisher_P);
				} else {
					throw new RuntimeException("p value is not in [0,1]!");
				}

			}
			return true;
		} else {
			return false;
		}
	}

	private List<GroupResult> writeGroupResult(List<RefCpG> refCpGList, ASMGraph graph) throws IOException {
		BufferedWriter groupResultWriter = new BufferedWriter(
				new FileWriter(inputFile.getAbsolutePath() + ".detected"));

		groupResultWriter.write(inputFile.getName() + "\n");
		groupResultWriter.write(String.format("tied weight counter:%d\n", graph.getTieWeightCounter()));
		groupResultWriter.write(String.format("tied id counter:%d\n", graph.getTieIdCountCounter()));

		printRefCpGGroupCoverage(refCpGList, graph, groupResultWriter);

		List<GroupResult> groupResultList = writeAlignedResult(refCpGList, graph, groupResultWriter);
		writePvalues(refCpGList, groupResultWriter);

		writeGroupReadList(groupResultWriter, groupResultList);

		groupResultWriter.close();

		return groupResultList;
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

	private void printRefCpGGroupCoverage(List<RefCpG> refCpGList, ASMGraph graph,
	                                      BufferedWriter groupResultWriter) throws IOException {
		// print refCpG's group coverage
		for (RefCpG refCpG : refCpGList) {
			if (graph.getClusterRefCpGMap().containsKey(refCpG.getPos())) {
				groupResultWriter.write(String.format("pos: %d\tcount: %d\n", refCpG.getPos(),
						graph.getClusterRefCpGMap().get(refCpG.getPos()).getClusterCount()));
			}
		}
	}

	private List<GroupResult> writeAlignedResult(List<RefCpG> refCpGList, ASMGraph graph,
	                                             BufferedWriter groupResultWriter) throws IOException {
		List<GroupResult> groupResultList = graph.getClusterResult()
				.values()
				.stream()
				.map(vertex -> new GroupResult(new ArrayList<>(vertex.getRefCpGMap().values()),
						vertex.getMappedReadList(), vertex.getMECScore()))
				.collect(Collectors.toList());

		// start to write aligned result
		groupResultWriter.write("Aligned result:\n");
		groupResultList.sort(GroupResult::compareTo);

		for (int i = 0; i < groupResultList.size(); i++) {
			groupResultWriter.write(i + ":\t");
			int count = 0;
			double methylLevel = 0;
			for (RefCpG refCpG : refCpGList) {
				Map<Integer, RefCpG> refCpGMap = groupResultList.get(i)
						.getRefCpGList()
						.stream()
						.collect(Collectors.toMap(RefCpG::getPos, r -> r));
				if (refCpGMap.containsKey(refCpG.getPos())) {
					count++;
					methylLevel += refCpGMap.get(refCpG.getPos()).getMethylLevel();
					groupResultWriter.write(Strings.padEnd(
							String.format("%.2f(%d)", refCpGMap.get(refCpG.getPos()).getMethylLevel(),
									refCpGMap.get(refCpG.getPos()).getCoveredCount()), ALIGN_COL_SIZE, ' '));
				} else {
					// fill the gap
					groupResultWriter.write(Strings.repeat(" ", ALIGN_COL_SIZE));
				}
			}
			groupResultList.get(i).setAvgMethylLevel(methylLevel / count);
			groupResultWriter.write("\n");
		}
		return groupResultList;
	}

	private void writePvalues(List<RefCpG> refCpGList, BufferedWriter groupResultWriter) throws IOException {
		groupResultWriter.write("P:\t");
		for (RefCpG refCpG : refCpGList) {
			if (refCpG.getP_value() != -1) {
				groupResultWriter.write(
						Strings.padEnd(String.format("%.4e", refCpG.getP_value()), ALIGN_COL_SIZE, ' '));
			} else {
				groupResultWriter.write(Strings.repeat(" ", ALIGN_COL_SIZE));

			}
		}
		groupResultWriter.write("\n");
	}

	private void writeGroupReadList(BufferedWriter groupResultWriter,
	                                List<GroupResult> groupResultList) throws IOException {
		groupResultWriter.write("Group's read list:\n");
		for (int i = 0; i < groupResultList.size(); i++) {
			groupResultWriter.write("Group:" + i + "\t" + "size:" + groupResultList.get(i).getMappedReadList().size() +
					"\tMEC:" +
					groupResultList.get(i).getMec() + "\n");
			for (MappedRead mappedRead : groupResultList.get(i).getMappedReadList()) {
				groupResultWriter.write(mappedRead.getId() + ",");
			}
			groupResultWriter.write("\n");
		}
	}

	private double EPCaluculation_average_combine(Vertex minorityCluster, Vertex majorityCluster,
	                                              List<RefCpG> twoClusterRefCpGList) {
		double p = 0;
		for (RefCpG refCpG : twoClusterRefCpGList) {
			double majorMethylLevel = majorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylLevel();
			double comb = CombinatoricsUtils.binomialCoefficientDouble(
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

	private double calcRegionP_StoufferComb(List<RefCpG> twoClusterRefCpGList) {
		NormalDistribution stdNorm = new NormalDistribution(0, 1);
		double z = twoClusterRefCpGList.stream()
				.mapToDouble(refCpG -> stdNorm.inverseCumulativeProbability(1 - refCpG.getP_value()))
				.sum() / Math.sqrt(twoClusterRefCpGList.size());
		return 1 - stdNorm.cumulativeProbability(z);
	}

	private double calcRegionP_FisherComb(List<RefCpG> twoClusterRefCpGList) {
		double regionP = -2 * twoClusterRefCpGList.stream().mapToDouble(refCpG -> Math.log(refCpG.getP_value())).sum();
		ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(2 * twoClusterRefCpGList.size());
		return 1 - chiSquaredDistribution.cumulativeProbability(regionP);
	}

	private double calcRegionP_SidakComb(List<RefCpG> twoClusterRefCpGList) {
		double minP = Double.MAX_VALUE;
		for (RefCpG refCpG : twoClusterRefCpGList) {
			if (minP > refCpG.getP_value()) {
				minP = refCpG.getP_value();
			}
		}
		return 1 - Math.pow(1 - minP, twoClusterRefCpGList.size());
	}
}
