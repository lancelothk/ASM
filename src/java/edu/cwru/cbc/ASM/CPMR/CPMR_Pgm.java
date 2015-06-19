package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.genomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * main entry of program.
 */
public class CPMR_Pgm {
	public static final int INIT_POS = 0;

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption("r", true, "Reference File");
		options.addOption("m", true, "MappedRead File");
		options.addOption("o", true, "Output Path");
		options.addOption("p", true, "Partial methylation threshold");
		options.addOption("mcc", true, "Minimum adjacent CpG coverage");
		options.addOption("mic", true, "Minimum interval CpG number");
		options.addOption("mir", true, "Minimum interval read number");

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceGenomeFileName = cmd.getOptionValue("r");
		String mappedReadFileName = cmd.getOptionValue("m");
		String outputPath = cmd.getOptionValue("o");
		double partial_methyl_threshold = Double.valueOf(cmd.getOptionValue("p"));
		int min_cpg_coverage = Integer.valueOf(cmd.getOptionValue("mcc"));
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		int min_interval_reads = Integer.valueOf(cmd.getOptionValue("mir"));

		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = IOUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		// load mapped reads
		start = System.currentTimeMillis();
		List<MappedRead> mappedReadList = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
				new MappedReadLineProcessor());
		System.out.println(
				"load mappedReadLinkedHashMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");
		Map<Integer, RefCpG> refMap = refCpGList.stream().collect(Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refMap));

		InputReadsSummary inputReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		mappedReadList.forEach(inputReadsSummary::addMappedRead);

		InputReadsSummary cpgReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		mappedReadList.forEach(mr -> {
			if (mr.getCpgList().size() > 0) {
				cpgReadsSummary.addMappedRead(mr);
			}
		});

		String summaryFileName = outputPath + "/CPMR.bed";
		String reportFileName = outputPath + "/CPMR.report";
		File intervalFolder = new File(outputPath);
		if (!intervalFolder.exists()) {
			if (!intervalFolder.mkdirs()) {
				throw new RuntimeException("failed to create folder:\t" + intervalFolder);
			}
		}

		CPMR cpmr = new CPMR(refCpGList, refChr, min_cpg_coverage, min_interval_cpg, min_interval_reads,
				partial_methyl_threshold);
		List<ImmutableGenomicInterval> immutableGenomicIntervals = cpmr.getGenomicIntervals();
		writeIntervals(outputPath, summaryFileName, immutableGenomicIntervals);
		InputReadsSummary intervalReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		immutableGenomicIntervals.forEach(i -> i.getMappedReadList().forEach(intervalReadsSummary::addMappedRead));
		List<RefCpG> refCpGCollection = immutableGenomicIntervals.stream().flatMap(
				i -> i.getRefCpGList().stream()).collect(Collectors.toList());
		writeReport(reportFileName,
				inputReadsSummary.getSummaryString(
						"Summary of raw input reads:\n") + cpgReadsSummary.getSummaryString(
						"\nSummary of reads with at least 1 CpG:\n") + intervalReadsSummary.getSummaryString(
						"\nSummary of reads in interval:\n") + InputReadsSummary.getCpGCoverageSummary(
						refCpGCollection), refCpGList.size(), cpmr.getRawIntervalCount(),
				immutableGenomicIntervals.size());
	}

	/**
	 * write single interval output in mapped read format like the original data *
	 */
	private static void writeMappedReadInInterval(String intervalFolderName,
	                                              ImmutableGenomicInterval immutableGenomicInterval) throws
			IOException {
		BufferedWriter mappedReadWriter = new BufferedWriter(
				new FileWriter(String.format("%s/%s-%d-%d%s", intervalFolderName, immutableGenomicInterval.getChr(),
						immutableGenomicInterval.getStart(), immutableGenomicInterval.getEnd(),
						Constant.MAPPEDREADS_EXTENSION)));
		mappedReadWriter.write(String.format("ref:\t%s\n",
				immutableGenomicInterval.getRefString().substring(immutableGenomicInterval.getStart() - INIT_POS,
						immutableGenomicInterval.getEnd() + 1 - INIT_POS)));
		for (MappedRead mappedRead : immutableGenomicInterval.getMappedReadList()) {
			mappedReadWriter.write(
					mappedRead.toMappedReadString(immutableGenomicInterval.getStart(),
							immutableGenomicInterval.getEnd()));
		}
		mappedReadWriter.close();
	}

	private static void writeIntervals(String intervalFolderName, String summaryFileName,
	                                   List<ImmutableGenomicInterval> immutableGenomicIntervalList) throws IOException {
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		intervalSummaryWriter.write("chr\tstart\tend\treadCount\tCpGCount\n");
		for (ImmutableGenomicInterval immutableGenomicInterval : immutableGenomicIntervalList) {
			try {
				intervalSummaryWriter.write(
						String.format("%s\t%d\t%d\t%d\t%d\n", immutableGenomicInterval.getChr(),
								immutableGenomicInterval.getStart(),
								immutableGenomicInterval.getEnd(), immutableGenomicInterval.getMappedReadList().size(),
								immutableGenomicInterval.getRefCpGList().size()));
				writeMappedReadInInterval(intervalFolderName, immutableGenomicInterval);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		intervalSummaryWriter.close();
	}

	private static void writeReport(String reportFileName, String reportString, int cpgCount, int rawIntervalCount,
	                                int outputIntervalCount) throws
			IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(reportFileName));
		writer.write(reportString);
		writer.write("Total #CpG in reference:\t" + cpgCount + "\n");
		writer.write("Raw Interval count:\t" + rawIntervalCount + "\n");
		writer.write("Output Interval count:\t" + outputIntervalCount + "\n");
		writer.close();
	}

}