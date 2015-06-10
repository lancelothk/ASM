package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.GenomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.IO.InputReadsSummary;
import edu.cwru.cbc.ASM.commons.IO.MappedReadLineProcessorWithFilter;
import edu.cwru.cbc.ASM.commons.Methylation.RefChr;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;

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
		options.addOption("mrc", true, "Minimum read CpG number");
		options.addOption("mic", true, "Minimum interval CpG number");
		options.addOption("mir", true, "Minimum interval read number");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceGenomeFileName = cmd.getOptionValue("r");
		String mappedReadFileName = cmd.getOptionValue("m");
		String outputPath = cmd.getOptionValue("o");
		double partial_methyl_threshold = Double.valueOf(cmd.getOptionValue("p"));
		int min_cpg_coverage = Integer.valueOf(cmd.getOptionValue("mcc"));
		int min_read_cpg = Integer.valueOf(cmd.getOptionValue("mrc"));
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		int min_interval_reads = Integer.valueOf(cmd.getOptionValue("mir"));

		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = CommonsUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		// load mapped reads
		start = System.currentTimeMillis();
		String reportString = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
				new MappedReadLineProcessorWithFilter(refCpGList, min_read_cpg, refChr.getRefString().length()));
		System.out.println("load mappedReadList complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		String summaryFileName = outputPath + "/CPMR.bed";
		String reportFileName = outputPath + "/CPMR.report";
		File intervalFolder = new File(outputPath);
		if (!intervalFolder.exists()) {
			if (!intervalFolder.mkdirs()) {
				throw new RuntimeException("failed to create folder:\t" + intervalFolder);
			}
		}

		CPMR cpmr = new CPMR(refCpGList, min_cpg_coverage, min_read_cpg,
				min_interval_cpg, min_interval_reads, partial_methyl_threshold);
		List<List<RefCpG>> cpgIntervalList = cpmr.getIntervals();
		List<ImmutableGenomicInterval> immutableGenomicIntervals = cpmr.getGenomicIntervals(refChr, cpgIntervalList);
		writeIntervals(outputPath, summaryFileName, immutableGenomicIntervals);
		InputReadsSummary intervalReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		immutableGenomicIntervals.forEach(i -> i.getMappedReadList().forEach(intervalReadsSummary::addMappedRead));
		List<RefCpG> refCpGCollection = immutableGenomicIntervals.stream().flatMap(
				i -> i.getRefCpGList().stream()).collect(Collectors.toList());
		writeReport(reportFileName,
				reportString + intervalReadsSummary.getSummaryString(
						"\nSummary of reads in interval:\n") + InputReadsSummary.getCpGCoverageSummary(
						refCpGCollection), cpgIntervalList.size(), immutableGenomicIntervals.size());
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
					mappedRead.outputString(immutableGenomicInterval.getStart(), immutableGenomicInterval.getEnd()));
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

	private static void writeReport(String reportFileName, String reportString, int rawIntervalCount,
	                                int outputIntervalCount) throws
			IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(reportFileName));
		writer.write(reportString);
		writer.write("Raw Interval count:\t" + rawIntervalCount + "\n");
		writer.write("Output Interval count:\t" + outputIntervalCount + "\n");
		writer.close();
	}

}