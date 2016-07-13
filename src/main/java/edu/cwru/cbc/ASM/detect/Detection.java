package edu.cwru.cbc.ASM.detect;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import edu.cwru.cbc.ASM.detect.dataType.ASMGraph;
import edu.cwru.cbc.ASM.detect.dataType.ClusterRefCpG;
import edu.cwru.cbc.ASM.detect.dataType.GroupResult;
import edu.cwru.cbc.ASM.detect.dataType.Vertex;
import edu.cwru.cbc.ASM.tools.ReadsVisualizationPgm;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;

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
public class Detection implements Callable<String> {
	private static final int ALIGN_COL_SIZE = 12;
	private File inputFile;
	private int startPos;
	private String chr;
	private int endPos;
	private int min_interval_cpg;
	private int permTime;
	private String outputPath;

	/**
	 * Detection constructor.
	 *  @param inputFile        A file contains reads in input region or a folder contains all input region files. File name should be in format:chr-start-end
	 * @param outputPath Path for output files
	 * @param min_interval_cpg Minimum number of CpGs in the interval. If under this threshold(too small), won't compute the error probability.
	 * @param permTime         Time of random permutation
	 */
	public Detection(File inputFile, String outputPath, int min_interval_cpg, int permTime) {
		this.inputFile = inputFile;
		this.min_interval_cpg = min_interval_cpg;
		this.permTime = permTime;
		this.outputPath = outputPath;
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


	public String call() throws IOException {
		String reference = IOUtils.readRefFromIntervalReadsFile(inputFile);
		extractIntervalPosition(inputFile);

		List<MappedRead> mappedReadList = readMappedReads(reference);
		ASMGraph partitionedGraph = executePartition(mappedReadList);

		Collection<Vertex> clusters = partitionedGraph.getClusterResult().values();
		alignPartitionedReads(clusters, reference);
		return "";
//
//
//		List<RefCpG> twoClusterRefCpGList = getTwoClusterRefCpGs(partitionedGraph);
//		if (clusters.size() > 2) {
//			// more than 2 clusters
//			throw new RuntimeException("more than two clusters!");
//		}
//		if (clusters.size() < 2) {
//			// less than 1 clusters
//			return "";
//		}
//		if (twoClusterRefCpGList.stream().count() < min_interval_cpg) {
//			// interval contain less #cpg than min_interval_cpg
//			return "";
//		}
//
//		// do tests
//		List<Double> pValueList = DetectionUtils.fisherExactTest(twoClusterRefCpGList, clusters);
//		double combSum = DetectionUtils.calcRegionP_FisherCombSum(pValueList);
//		double regionP = 1 - new ChiSquaredDistribution(2 * twoClusterRefCpGList.size()).cumulativeProbability(combSum);
//		double clusterIndex = DetectionUtils.calcClusterIndex(twoClusterRefCpGList, clusters);
//
//		// index of random data has larger score
//		int randIndex = getRandomIndex(combSum, reference);
//
//		if (randIndex >= 0) {
//			return "";
//		}
//
//		List<GroupResult> groupResultList = writePartitionResult(clusters, reference, twoClusterRefCpGList);
//
//		return IntervalDetectionSummaryFormatter.formatSummaryString(chr, startPos, endPos, endPos - startPos + 1,
//				startPos + "-" + endPos,
//				partitionedGraph.getOriginalEdgeCount(), mappedReadList.size(), twoClusterRefCpGList.size(),
//				clusters.size(), partitionedGraph.getCpGSum(), partitionedGraph.getMECSum(),
//				partitionedGraph.getNormMECSum(), regionP, combSum, clusterIndex,
//				Iterables.get(clusters, 0).getMappedReadList().size(),
//				Iterables.get(clusters, 1).getMappedReadList().size(),
//				groupResultList.get(0).getAvgMethylLevel(),
//				groupResultList.get(1).getAvgMethylLevel());
	}

	private int getRandomIndex(double combSum, String reference) throws IOException {
		List<MappedRead> mappedReadList = readMappedReads(reference);
		for (int i = 0; i <= permTime; i++) {
			ASMGraph randGraph = executePartition(randomizeMethylStatus(mappedReadList));
			List<RefCpG> twoClusterRefCpGList = getTwoClusterRefCpGs(randGraph);
			List<Double> pValueList = DetectionUtils.fisherExactTest(twoClusterRefCpGList,
					randGraph.getClusterResult().values());
			double randCombSum = DetectionUtils.calcRegionP_FisherCombSum(pValueList);
			if (combSum < randCombSum) {
				return i;
			}
		}
		return -1;
	}

	private ASMGraph executePartition(List<MappedRead> mappedReadList) {
		ASMGraph randGraph = new ASMGraph(mappedReadList);
		randGraph.cluster();
		return randGraph;
	}

	private List<MappedRead> readMappedReads(String reference) throws IOException {
		// load input
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		// filter out reads which covers no CpG sites
		return Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor(mr -> mr.generateCpGsInRead(refMap) > 0));
	}

	private List<RefCpG> getTwoClusterRefCpGs(ASMGraph graph) {
		return graph.getClusterRefCpGMap()
				.values()
				.stream()
				.filter(clusterRefCpG -> clusterRefCpG.getClusterCount() == 2)
				.map(ClusterRefCpG::getRefCpG)
				.sorted()
				.collect(Collectors.toList());
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

	public List<GroupResult> writePartitionResult(Collection<Vertex> clusters, String reference,
	                                              List<RefCpG> refCpGList) throws IOException {
		// *.detected
		List<GroupResult> groupResultList = writeDetectedResults(clusters, refCpGList);
		// *.groups.aligned
		alignPartitionedReads(clusters, reference);
		return groupResultList;
	}

	private List<GroupResult> writeDetectedResults(Collection<Vertex> clusters, List<RefCpG> refCpGList) throws
			IOException {
		BufferedWriter groupResultWriter = new BufferedWriter(
				new FileWriter(outputPath + "/" + inputFile.getName() + ".detected"));
		List<GroupResult> groupResultList = writeAlignedResult(refCpGList, groupResultWriter, clusters);
		writePvalues(refCpGList, groupResultWriter);
		writeGroupReadList(groupResultWriter, groupResultList);
		groupResultWriter.close();
		return groupResultList;
	}

	private void alignPartitionedReads(Collection<Vertex> clusterResults, String reference) throws IOException {
		List<List<MappedRead>> readGroups = clusterResults.stream()
				.map(Vertex::getMappedReadList)
				.collect(Collectors.toList());
		ReadsVisualizationPgm.writeAlignedReadsIntoGroups(readGroups, reference,
				outputPath + "/" + inputFile.getName() + ".groups.aligned");
	}

	private List<GroupResult> writeAlignedResult(List<RefCpG> refCpGList, BufferedWriter groupResultWriter,
	                                             Collection<Vertex> clusterResults) throws IOException {
		List<GroupResult> groupResultList = clusterResults.stream()
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
		for (int i = 0; i < groupResultList.size(); i++) {
			groupResultWriter.write("Group:" + i + "\t" + "size:" + groupResultList.get(i).getMappedReadList().size() +
					"\tMEC:" + groupResultList.get(i).getMec() + "\tavgMethyl:" + groupResultList.get(i)
					.getAvgMethylLevel() + "\n");
		}
	}
}
