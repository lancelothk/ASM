package edu.cwru.cbc.ASM.detect.WithMappedRead;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.align.AlignReads;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import edu.cwru.cbc.ASM.detect.Detection;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.ASMGraph;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.ClusterRefCpG;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.GroupResult;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.Vertex;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.util.FastMath;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;

/**
 * Created by kehu on 11/12/14.
 * ASM Detection with whole read info.
 */
public class DetectionWithMappedRead extends Detection {
	private static final Logger logger = Logger.getLogger(DetectionWithMappedRead.class.getName());
	private static final String EXPERIMENT_NAME = "2group";
	private static int countOfCpGLessThanFive = 0;

	public DetectionWithMappedRead() {
	}

	public static void main(
			String[] args) throws ExecutionException, InterruptedException, IOException, InstantiationException, IllegalAccessException {
		long start = System.currentTimeMillis();
		String cellLine = "i90";
		String replicate = "r1";
		String name = "chr20";
		String homeDirectory = System.getProperty("user.home");

		final int MIN_CONT_COVERAGE = 4;
		final int MIN_INTERVAL_READS = 10;
		final int MIN_INTERVAL_CPG = 5;

		String fileName = "";
		String inputName = String.format("%s/experiments/ASM/result_%s_%s/intervals_%s_%d_%d_%d/%s", homeDirectory,
										 cellLine, replicate, name, MIN_CONT_COVERAGE, MIN_INTERVAL_CPG,
										 MIN_INTERVAL_READS, fileName);
		String summaryFileName = String.format(
				"%1$s/experiments/ASM/result_%2$s_%3$s/%2$s_%3$s_%4$s_ASM_summary_%5$s_%6$s_%7$d_%8$d_%9$d",
				homeDirectory, cellLine, replicate, name, fileName, EXPERIMENT_NAME, MIN_CONT_COVERAGE,
				MIN_INTERVAL_CPG, MIN_INTERVAL_READS);

		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		summaryWriter.write(
				"chr\tstartPos\tendPos\tlength\tvertex number\tedge number\treadCount\tCpGCount\tGroupCount\t" +
						"avgGroupPerCpG\tMECSum\tCpGSum\tNormMEC\tmwwScore\tmwwTest\terrorProbability\tgroupSizes\n");

		int threadNumber = 6;

		List<String> resultList = new DetectionWithMappedRead().execute(inputName, threadNumber);

		for (String result : resultList) {
			summaryWriter.write(result);
		}

		summaryWriter.close();
		logger.info("count of intervals with overlappedCpG < 5\t" + countOfCpGLessThanFive);
		logger.info(System.currentTimeMillis() - start + "ms");
	}

	public String call() throws Exception {
		extractIntervalPosition(inputFile);

		// load input
		String reference = Utils.readRefFromIntervalFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		// filter out reads which only cover 1 or no CpG sites
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8).readLines(
				new MappedReadLineProcessor(refCpGList));

		// align read
		AlignReads.alignReads(mappedReadList, reference, inputFile.getAbsolutePath() + ".aligned");

		// construct graph
		ASMGraph graph = new ASMGraph(mappedReadList);

		// clustering
		graph.cluster();

		double avgGroupCpGCoverage = 0;//writeGroupResult(refCpGList, graph, inputFile);

		return buildSummary(endPos - startPos + 1, refCpGList, mappedReadList, graph, avgGroupCpGCoverage);
	}

	private double calcErrorProbability(List<RefCpG> refCpGList, Collection<Vertex> clusterResult,
										Map<Integer, ClusterRefCpG> clusterRefCpGMap) {
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

				List<RefCpG> twoClusterRefCpGList = new ArrayList<>();
				for (RefCpG refCpG : refCpGList) {
					if (clusterRefCpGMap.containsKey(refCpG.getPos())) {
						if (clusterRefCpGMap.get(refCpG.getPos()).getClusterCount() == 2) {
							twoClusterRefCpGList.add(refCpG);
						}
					}
				}

				// overlapped CpG < 5
				if (twoClusterRefCpGList.size() < 5) {
					p = -3;
				} else {
					p = 0;
					for (RefCpG refCpG : twoClusterRefCpGList) {
						double majorMethylLevel = majorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylLevel();
						p += FastMath.pow(majorMethylLevel,
										  minorityCluster.getRefCpGMap().get(refCpG.getPos()).getMethylCount()) *
								FastMath.pow(1 - majorMethylLevel,
											 minorityCluster.getRefCpGMap().get(refCpG.getPos()).getNonMethylCount());
					}
					p /= twoClusterRefCpGList.size();
				}
				break;
			default:
				p = -2; // more than 2 clusters
		}
		return p;
	}

	private MWW calculateMWW(List<RefCpG> refCpGList, ASMGraph graph) {
		MWW MWW = new MWW();
		if (graph.getClusterResult().size() == 1) {
			MWW.mwwScore = -1; // 1 cluster interval
			MWW.mwwPvalue = -1;
		} else { // cluster size == 2 or more
			List<RefCpG> twoClusterRefCpGList = new ArrayList<>();
			for (RefCpG refCpG : refCpGList) {
				if (graph.getClusterRefCpGMap().containsKey(refCpG.getPos())) {
					if (graph.getClusterRefCpGMap().get(refCpG.getPos()).getClusterCount() == 2) {
						twoClusterRefCpGList.add(refCpG);
					}
				}
			}
			// overlapped CpG < 5
			if (twoClusterRefCpGList.size() < 5) {
				countOfCpGLessThanFive++;
			}
			MannWhitneyUTest test = new MannWhitneyUTest();
			double[][] mwwMatrix = new double[2][twoClusterRefCpGList.size()];

			for (int i = 0; i < twoClusterRefCpGList.size(); i++) {
				assert twoClusterRefCpGList.get(i).getCpGCoverage() == 2;
				int j = 0;
				for (Vertex vertex : graph.getClusterResult().values()) {
					if (vertex.getRefCpGMap().containsKey(twoClusterRefCpGList.get(i).getPos())) {
						mwwMatrix[j][i] = vertex.getRefCpGMap().get(
								twoClusterRefCpGList.get(i).getPos()).getMethylLevel();
						j++;
					}
				}
			}
			MWW.mwwScore = test.mannWhitneyU(mwwMatrix[0], mwwMatrix[1]);
			MWW.mwwPvalue = test.mannWhitneyUTest(mwwMatrix[0], mwwMatrix[1]);
		}
		return MWW;

	}

	private String buildSummary(int length, List<RefCpG> refCpGList, List<MappedRead> mappedReadList, ASMGraph graph,
								double avgGroupCpGCoverage) {
		MWW MWW = calculateMWW(refCpGList, graph);

		StringBuilder sb = new StringBuilder(
				String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\t%d\t%f\t%f\t%f\t%f\t", chr, startPos, endPos,
							  length, graph.getOriginalVertexCount(), graph.getOriginalEdgeCount(),
							  mappedReadList.size(), refCpGList.size(), graph.getClusterResult().size(),
							  avgGroupCpGCoverage, graph.getMECSum(), graph.getCpGSum(), graph.getNormMECSum(),
							  MWW.mwwScore, MWW.mwwPvalue,
							  calcErrorProbability(refCpGList, graph.getClusterResult().values(),
												   graph.getClusterRefCpGMap())));
		for (Vertex vertex : graph.getClusterResult().values()) {
			sb.append(vertex.getMappedReadList().size()).append(",");
		}
		sb.append("\n");
		return sb.toString();
	}

	private double writeGroupResult(List<RefCpG> refCpGList, ASMGraph graph, File inputFile) throws IOException {
		BufferedWriter groupResultWriter = new BufferedWriter(
				new FileWriter(inputFile.getAbsolutePath() + "." + EXPERIMENT_NAME));

		groupResultWriter.write(inputFile.getName() + "\n");
		groupResultWriter.write(String.format("tied weight counter:%d\n", graph.getTieWeightCounter()));
		groupResultWriter.write(String.format("tied id counter:%d\n", graph.getTieIdCountCounter()));

		// assign index of cpg
		for (int i = 0; i < refCpGList.size(); i++) {
			refCpGList.get(i).assignIndex(i);
		}

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

		List<GroupResult> groupResultList = graph.getClusterResult().values().stream().map(
				vertex -> new GroupResult(new ArrayList<>(vertex.getRefCpGMap().values()), vertex.getMappedReadList(),
										  vertex.getMECScore())).collect(Collectors.toList());

		// start to write aligned result
		groupResultWriter.write("Aligned result:\n");
		groupResultList.sort(GroupResult::compareTo);

		for (int i = 0; i < groupResultList.size(); i++) {
			groupResultWriter.write(i + "\t");
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

		groupResultWriter.close();
		return sum / graph.getClusterRefCpGMap().size();
	}

	@Override
	protected Detection constructNewInstance() {
		return new DetectionWithMappedRead();
	}

	private class MWW {
		public double mwwScore = -1010, mwwPvalue = -1010;
	}
}
