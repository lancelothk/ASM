package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.CpG.CpG;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;
import edu.cwru.cbc.ASM.commons.Read.MappedReadLineProcessorWithSummary;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * CPMR stands for Continuous covered partial methylation region.
 * Main entry of the module.
 */
public class CPMR {
	public static final int INIT_POS = 0;
	private static int outputIntervalCount = 0;

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption("r", true, "Reference File");
		options.addOption("m", true, "MappedRead File");
		options.addOption("o", true, "Output Path");
		options.addOption("mcc", true, "Minimum CpG coverage");
		options.addOption("mrc", true, "Minimum read CpG number");
		options.addOption("mic", true, "Minimum interval CpG number");
		options.addOption("mir", true, "Minimum interval read number");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceGenomeFileName = cmd.getOptionValue("r");
		String mappedReadFileName = cmd.getOptionValue("m");
		String outputPath = cmd.getOptionValue("o");
		int min_cpg_coverage = Integer.valueOf(cmd.getOptionValue("mcc"));
		int min_read_cpg = Integer.valueOf(cmd.getOptionValue("mrc"));
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		int min_interval_reads = Integer.valueOf(cmd.getOptionValue("mir"));

		splitEpigenome(referenceGenomeFileName, mappedReadFileName, outputPath, min_cpg_coverage, min_read_cpg,
				min_interval_cpg, min_interval_reads);
	}

	private static void splitEpigenome(String referenceGenomeFileName, String mappedReadFileName, String outputPath,
									   int min_cpg_coverage, int min_read_cpg, int min_interval_cpg,
									   int min_interval_reads) throws IOException {
		long start = System.currentTimeMillis();
		outputPath += String.format("/experiment_%d_%d_%d_%d/", min_cpg_coverage, min_read_cpg, min_interval_cpg,
				min_interval_reads);
		String intervalFolderName = outputPath + "intervals";
		String summaryFileName = outputPath + "CPMR_summary";
		String reportFileName = outputPath + "CPMR_report";
		File intervalFolder = new File(intervalFolderName);
		if (!intervalFolder.exists()) {
			//noinspection ResultOfMethodCallIgnored
			intervalFolder.mkdirs();
		}

		RefChr refChr = CommonsUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		refCpGList.sort((RefCpG c1, RefCpG c2) -> c1.getPos() - c2.getPos());
		System.out.println("cpg sorting complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		// assign cpg id after sorted
		for (int i = 0; i < refCpGList.size(); i++) {
			refCpGList.get(i).assignIndex(i);
		}
		System.out.println("cpg id assignment complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
				new MappedReadLineProcessorWithSummary(refCpGList, min_read_cpg, refChr.getRefString().length(),
						reportFileName));
		System.out.println("load mappedReadList complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		List<List<RefCpG>> cpgSiteIntervalList = getIntervals(refCpGList, min_cpg_coverage);

		writeIntervals(intervalFolderName, summaryFileName, refChr, cpgSiteIntervalList, min_cpg_coverage,
				min_interval_cpg, min_interval_reads, min_read_cpg);
		writeReport(reportFileName, cpgSiteIntervalList);
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");
	}

	private static List<List<RefCpG>> getIntervals(List<RefCpG> refCpGList, int min_cpg_coverage) {
		// split regions by continuous CpG coverage
		List<List<RefCpG>> cpgSiteIntervalList = new ArrayList<>();
		boolean cont = true;
		List<RefCpG> resultRefCpGList = new ArrayList<>();
		for (int i = 0; i < refCpGList.size() - 1; i++) {
			RefCpG curr = refCpGList.get(i);
			if (curr.getCpGCoverage() >= min_cpg_coverage && curr.hasPartialMethyl()) {
				if (!cont && resultRefCpGList.size() > 0) {
					cpgSiteIntervalList.add(resultRefCpGList);
					resultRefCpGList = new ArrayList<>();
				}
				resultRefCpGList.add(curr);
				RefCpG next = refCpGList.get(i + 1);
				cont = next.getCpGCoverage() >= min_cpg_coverage && next.hasPartialMethyl() && curr.hasCommonRead(next);
			}
		}
		if (resultRefCpGList.size() > 0) {
			if (cont) {
				resultRefCpGList.add(resultRefCpGList.get(resultRefCpGList.size() - 1));
			}
			cpgSiteIntervalList.add(resultRefCpGList);
		}
		return cpgSiteIntervalList;
	}

	private static void writeIntervals(String intervalFolderName, String summaryFileName, RefChr refChr,
									   List<List<RefCpG>> cpgSiteIntervalList, int min_cpg_coverage,
									   int min_interval_cpg, int min_interval_reads,
									   int min_read_cpg) throws IOException {
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		intervalSummaryWriter.write("chr\tlength\treadCount\tCpGCount\tstartCpG\tendCpG\n");
		cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = getReadsFromCpGs(min_cpg_coverage, list);
			// only pass high quality result for next step.
			if (mappedReadSet.size() >= min_interval_reads && list.size() >= min_interval_cpg) {
				List<MappedRead> mappedReadList = mappedReadSet.stream()
						.sorted((m1, m2) -> m1.getId().compareTo(m2.getId()))
						.sorted((m1, m2) -> m1.getStart() - m2.getStart())
						.collect(Collectors.toList());
				outputIntervalCount++;
				int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				@SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
						int startPos = startCpGPos;
				int endPos = endCpGPos + 1;
				try {
					intervalSummaryWriter.write(
							String.format("%s\t%d\t%d\t%d\t%d\t%d\n", refChr.getChr(), endPos - startPos + 1,
									mappedReadList.size(), list.size(), startCpGPos, endCpGPos));
					writeMappedReadInInterval(intervalFolderName, refChr, startPos, endPos, mappedReadList,
							min_read_cpg);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		});
		intervalSummaryWriter.close();
	}

	private static void writeReport(String reportFileName, List<List<RefCpG>> cpgSiteIntervalList) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(reportFileName, true));
		writer.write("Raw Interval count:\t" + cpgSiteIntervalList.size() + "\n");
		writer.write("Output Interval count:\t" + outputIntervalCount + "\n");
		writer.close();
	}

	private static Set<MappedRead> getReadsFromCpGs(int min_cpg_coverage, List<RefCpG> list) {
		Set<MappedRead> mappedReadSet = new HashSet<>();
		list.forEach((refCpG) -> {
			refCpG.getCpGList().forEach((CpG cpg) -> mappedReadSet.add(cpg.getMappedRead()));
			assert refCpG.getCpGList().size() >= min_cpg_coverage; // make sure coverage is correct
		});
		return mappedReadSet;
	}

	/**
	 * write single interval output in mapped read format like the original data *
	 */
	public static void writeMappedReadInInterval(String intervalFolderName, RefChr refChr, int startPos, int endPos,
												 Collection<MappedRead> mappedReadSet,
												 int min_read_cpg) throws IOException {
		BufferedWriter mappedReadWriter = new BufferedWriter(
				new FileWriter(String.format("%s/%s-%d-%d", intervalFolderName, refChr.getChr(), startPos, endPos)));
		mappedReadWriter.write(String.format("ref:\t%s\n",
				refChr.getRefString().substring(startPos - INIT_POS, endPos + 1 - INIT_POS)));
		for (MappedRead mappedRead : mappedReadSet) {
			// The mapped read to write should contain at least min_read_cpg cpgs.
			assert mappedRead.getCpgList().size() < min_read_cpg;
			mappedReadWriter.write(mappedRead.outputString(startPos, endPos));
		}
		mappedReadWriter.close();
	}
}