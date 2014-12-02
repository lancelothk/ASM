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
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.GroupResult;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DataType.Vertex;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;

/**
 * Created by kehu on 11/12/14.
 * ASM Detection with whole read info.
 */
public class DetectionWithMappedRead extends Detection {
	private static final int MIN_READ_CPG = 2;

	public DetectionWithMappedRead() {
	}

	public static void main(
			String[] args) throws ExecutionException, InterruptedException, IOException, InstantiationException, IllegalAccessException {
		long start = System.currentTimeMillis();
		// test reads setting
//        String pathName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/chrTest2-1-6";
//        String summaryFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.summary";
//        String groupResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.groupResult";

		String cellLine = "i90";
		String replicate = "r1";
		String name = "chr20";
		String homeDirectory = System.getProperty("user.home");

		String fileName = "";
//        String fileName = "";
		String inputName = String.format("%s/experiments/ASM/result_%s_%s/intervals_%s/%s", homeDirectory, cellLine,
										 replicate, name, fileName);
		String summaryFileName = String.format(
				"%1$s/experiments/ASM/result_%2$s_%3$s/%2$s_%3$s_%4$s_ASM_summary_%5$s_test", homeDirectory, cellLine,
				replicate, name, fileName);

		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		summaryWriter.write(
				"name\tlength\tvertex number\tedge number\treadCount\tCpGCount\tGroupCount\tavgGroupPerCpG\tMECSum\tCpGSum\tgroupSizes\n");

		int threadNumber = 6;

		List<String> resultList = new DetectionWithMappedRead().execute(inputName, threadNumber);

		for (String result : resultList) {
			summaryWriter.write(result);
		}

		summaryWriter.close();
		System.out.println(System.currentTimeMillis() - start + "ms");
	}

	public String call() throws Exception {
		String[] items = inputFile.getName().split("-");
		if (items.length != 3) {
			throw new RuntimeException("invalid input file name format!\t" + inputFile.getName());
		}
		int initPos = Integer.parseInt(items[1]);

		// load input
		String reference = readRef(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, initPos);
		// filter out reads which only cover 1 or no CpG sites
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8).readLines(
				new MappedReadLineProcessor(refCpGList, read -> read.getCpgList().size() >= MIN_READ_CPG));

		// align read
		AlignReads.alignReads(mappedReadList, reference, inputFile.getAbsolutePath() + ".aligned");

		// construct graph
		ASMGraph graph = new ASMGraph(mappedReadList);

		// clustering
		graph.getClusters();

		if (graph.getClusterResult().size() == 0) {
			System.out.printf("");
		}

		double avgGroupCpGCoverage = writeGroupResult(refCpGList, graph, inputFile);

		return buildSummary(Integer.parseInt(items[2]) - Integer.parseInt(items[1]) + 1, refCpGList, mappedReadList,
							graph, avgGroupCpGCoverage);
	}

	private String buildSummary(int length, List<RefCpG> refCpGList, List<MappedRead> mappedReadList, ASMGraph graph,
								double avgGroupCpGCoverage) {
		StringBuilder sb = new StringBuilder(
				String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\t%d\t", inputFile.getName(), length,
							  graph.getOriginalVertexCount(), graph.getOriginalEdgeCount(), mappedReadList.size(),
							  refCpGList.size(), graph.getClusterResult().size(), avgGroupCpGCoverage,
							  graph.getMECSum(), graph.getCpGSum()));
		for (Vertex vertex : graph.getClusterResult().values()) {
			sb.append(vertex.getMappedReadList().size()).append(",");
		}
		sb.append("\n");
		return sb.toString();
	}

	private double writeGroupResult(List<RefCpG> refCpGList, ASMGraph graph, File inputFile) throws IOException {
		BufferedWriter groupResultWriter = new BufferedWriter(
				new FileWriter(inputFile.getAbsolutePath() + ".new_group"));

		groupResultWriter.write(inputFile.getName() + "\n");
		groupResultWriter.write(String.format("tied weight counter:%d\n", graph.getTieWeightCounter()));
		groupResultWriter.write(String.format("tied methyl polarity:%d\n", graph.getTieMethylPolarity()));
		groupResultWriter.write(String.format("tied id counter:%d\n", graph.getTieIdCountCounter()));

		// assign index of cpg
		for (int i = 0; i < refCpGList.size(); i++) {
			refCpGList.get(i).assignIndex(i);
		}

		// print refCpG's group coverage
		double sum = 0;
		for (RefCpG refCpG : refCpGList) {
			if (graph.getCoveredCpGMap().containsKey(refCpG.getPos())) {
				groupResultWriter.write(String.format("pos: %d\tcount: %d\n", refCpG.getPos(),
													  graph.getCoveredCpGMap().get(refCpG.getPos())));
				sum += graph.getCoveredCpGMap().get(refCpG.getPos());
			}
		}

		List<GroupResult> groupResultList = graph.getClusterResult().values().stream().map(
				vertex -> new GroupResult(new ArrayList<>(vertex.getRefCpGMap().values()), vertex.getMappedReadList(),
										  vertex.getMECScore())).collect(Collectors.toList());

		// start to write aligned result
		groupResultWriter.write("Aligned result:\n");
		groupResultList.sort(GroupResult::compareTo);
		Map<Integer, Integer> posIndexMap = refCpGList.stream().collect(
				Collectors.toMap(RefCpG::getPos, RefCpG::getIndex));

		int startIndex = posIndexMap.get(groupResultList.get(0).getFirstPos());
		for (int i = 0; i < groupResultList.size(); i++) {
			groupResultWriter.write(i + "\t");
			int currIndex = posIndexMap.get(groupResultList.get(i).getFirstPos());
			// fill the gap
			for (int j = 0; j < currIndex - startIndex; j++) {
				groupResultWriter.write(Strings.repeat(".", 10));
			}
			for (RefCpG refCpG : groupResultList.get(i).getRefCpGList()) {
				groupResultWriter.write(
						Strings.padEnd(String.format("%.2f(%d)", refCpG.getMethylLevel(), refCpG.getCoveredCount()), 10,
									   ' '));
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
		return sum / graph.getCoveredCpGMap().size();
	}

	private String readRef(File inputFile) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFile));
		String line = bufferedReader.readLine();
		String[] items = line.split("\t");
		if (items.length != 2) {
			throw new RuntimeException(
					"invalid reference line in interval read file!\t" + line + "\t" + inputFile.getName());
		} else {
			return items[1].toUpperCase();
		}
	}

	@Override
	protected Detection constructNewInstance() {
		return new DetectionWithMappedRead();
	}
}
