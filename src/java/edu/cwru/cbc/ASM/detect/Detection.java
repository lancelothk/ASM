package edu.cwru.cbc.ASM.detect;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import edu.cwru.cbc.ASM.detect.dataType.*;
import edu.cwru.cbc.ASM.tools.ReadsVisualizationPgm;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by lancelothk on 6/8/15.
 * ASM Detection with whole read info.
 */
public class Detection implements Callable<IntervalDetectionSummary> {
	private static final int ALIGN_COL_SIZE = 12;
	private File inputFile;
	private String chr;
	private int startPos;
	private int endPos;
	private int min_interval_cpg;
	private int permTime;

	/**
	 * Detection constructor.
	 *
	 * @param inputFile        A file contains reads in input region or a folder contains all input region files. File name should be in format:chr-start-end
	 * @param min_interval_cpg Minimum number of CpGs in the interval. If under this threshold(too small), won't compute the error probability.
	 * @param permTime         Time of random permutation
	 */
	public Detection(File inputFile, int min_interval_cpg, int permTime) {
		this.inputFile = inputFile;
		this.min_interval_cpg = min_interval_cpg;
		this.permTime = permTime;
	}

	public IntervalDetectionSummary call() throws Exception {
		extractIntervalPosition(inputFile);

		// load input
		String reference = IOUtils.readRefFromIntervalReadsFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		// filter out reads which covers no CpG sites
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor(mr -> mr.generateCpGsInRead(refMap) > 0));

		// construct graph
		ASMGraph graph = new ASMGraph(mappedReadList);

		// clustering
		graph.cluster();

		List<RefCpG> twoClusterRefCpGList = getTwoClustersRefCpG(refCpGList, graph.getClusterRefCpGMap());
		Collection<Vertex> clusterResults = graph.getClusterResult().values();

		DetectionUtils.fisherTest(twoClusterRefCpGList, clusterResults);
		double regionP;
		double combSum = 0;
		if (clusterResults.size() > 2) {
			// more than 2 clusters
			throw new RuntimeException("more than two clusters!");
		} else if (clusterResults.size() < 2) {
			// less than 2 clusters
			regionP = 2;
		} else if (twoClusterRefCpGList.stream()
				.count() < min_interval_cpg) {
			// interval contain less #cpg than min_interval_cpg
			regionP = 3;
		} else {
			combSum = DetectionUtils.calcRegionP_FisherCombSum(twoClusterRefCpGList);
			ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(2 * twoClusterRefCpGList.size());
			regionP = 1 - chiSquaredDistribution.cumulativeProbability(combSum);
		}

		// calculate daviesBouldin index
		double clusterIndex = DetectionUtils.calcClusterIndex(twoClusterRefCpGList, graph.getClusterResult().values());

		List<GroupResult> groupResultList = writeGroupResult(inputFile.getAbsolutePath(), refCpGList, graph);
		if (graph.getClusterResult().values().size() != 2) {
			throw new RuntimeException("clusters number is not 2!");
		}
		List<List<MappedRead>> readGroups = graph.getClusterResult().values().stream().map(
				Vertex::getMappedReadList).collect(
				Collectors.toList());
		ReadsVisualizationPgm.writeAlignedReadsIntoGroups(readGroups, reference,
				inputFile.getAbsolutePath() + ".groups.aligned");

		// random group P value
		double mec = graph.getMECSum();
		double normMEC = graph.getNormMECSum();
		int minPCount = 0;
		if (regionP <= 1 && regionP >= 0) {
			minPCount = getRandScoreIndex(refCpGList, mappedReadList, combSum);
		}

		return new IntervalDetectionSummary(regionP, minPCount > 0, chr.replace("chr", ""), startPos, endPos,
				endPos - startPos + 1, startPos + "-" + endPos,
				graph.getOriginalEdgeCount(), mappedReadList.size(), refCpGList.size(),
				twoClusterRefCpGList.size(), graph.getClusterResult().size(), graph.getCpGSum(),
				mec, normMEC, regionP, minPCount, clusterIndex,
				Iterables.get(graph.getClusterResult().values(), 0).getMappedReadList().size(),
				Iterables.get(graph.getClusterResult().values(), 1).getMappedReadList().size(),
				groupResultList.get(0).getAvgMethylLevel(),
				groupResultList.get(1).getAvgMethylLevel(),
				"<label>");
	}

	private int getRandScoreIndex(List<RefCpG> refCpGList, List<MappedRead> mappedReadList, double combSum) {
		for (int i = 1; i <= permTime; i++) {
			ASMGraph randGraph = new ASMGraph(randomizeMethylStatus(mappedReadList));
			randGraph.cluster();
			List<RefCpG> twoClusterRefCpGList = getTwoClustersRefCpG(refCpGList, randGraph.getClusterRefCpGMap());
			DetectionUtils.fisherTest(twoClusterRefCpGList, randGraph.getClusterResult().values());
			double randCombSum = DetectionUtils.calcRegionP_FisherCombSum(twoClusterRefCpGList);
			if (combSum < randCombSum) {
				return i;
			}
		}
		return 0;
	}

	private List<MappedRead> randomizeMethylStatus(List<MappedRead> mappedReadList) {
		Random rand = new Random();
		for (MappedRead mappedRead : mappedReadList) {
			for (CpG cpG : mappedRead.getCpgList()) {
				cpG.setMethylStatus(rand.nextBoolean() ? MethylStatus.C : MethylStatus.T);
			}
		}
		return mappedReadList;
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
		return refCpGList.stream().filter(
				refCpG -> clusterRefCpGMap.get(refCpG.getPos())
						.getClusterCount() == 2 && clusterRefCpGMap.containsKey(
						refCpG.getPos())).collect(Collectors.toList());
	}

	private List<GroupResult> writeGroupResult(String inputFilePath, List<RefCpG> refCpGList, ASMGraph graph) throws
			IOException {
		BufferedWriter groupResultWriter = new BufferedWriter(
				new FileWriter(inputFilePath + ".detected"));

		groupResultWriter.write(inputFilePath + "\n");
		groupResultWriter.write(String.format("tied weight counter:%d\n", graph.getTieWeightCounter()));
		groupResultWriter.write(String.format("tied id counter:%d\n", graph.getTieIdCountCounter()));

		printRefCpGGroupCoverage(refCpGList, groupResultWriter, graph.getClusterRefCpGMap());

		List<GroupResult> groupResultList = writeAlignedResult(refCpGList, groupResultWriter, graph.getClusterResult()
				.values());
		writePvalues(refCpGList, groupResultWriter);

		writeGroupReadList(groupResultWriter, groupResultList);

		groupResultWriter.close();

		return groupResultList;
	}

	private void printRefCpGGroupCoverage(List<RefCpG> refCpGList, BufferedWriter groupResultWriter,
	                                      Map<Integer, ClusterRefCpG> clusterRefCpGMap) throws IOException {
		// print refCpG's group coverage
		for (RefCpG refCpG : refCpGList) {
			if (clusterRefCpGMap.containsKey(refCpG.getPos())) {
				groupResultWriter.write(String.format("pos: %d\tcount: %d\n", refCpG.getPos(),
						clusterRefCpGMap.get(refCpG.getPos()).getClusterCount()));
			}
		}
	}

	private List<GroupResult> writeAlignedResult(List<RefCpG> refCpGList,
	                                             BufferedWriter groupResultWriter,
	                                             Collection<Vertex> clusterResults) throws IOException {
		List<GroupResult> groupResultList = clusterResults
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
}
