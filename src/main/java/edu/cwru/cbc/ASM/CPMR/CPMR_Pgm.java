package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CMDHelper;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.MappedReadFileFormat;
import edu.cwru.cbc.ASM.commons.genomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * main entry of program.
 */
public class CPMR_Pgm {

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("r").hasArg().desc("Reference File (Required)").required().build());
		options.addOption(Option.builder("m").hasArg().desc("MappedRead File (Required)").required().build());
		options.addOption(Option.builder("o").hasArg().desc("Output Path (Required)").required().build());
		options.addOption(
				Option.builder("p").hasArg().desc("Partial methylation threshold (Required)").required().build());
		options.addOption(
				Option.builder("mcc").hasArg().desc("Minimum adjacent CpG coverage (Required)").required().build());
		options.addOption(
				Option.builder("mic").hasArg().desc("Minimum interval CpG number (Required)").required().build());
		options.addOption(
				Option.builder("mir").hasArg().desc("Minimum interval read number (Required)").required().build());
		options.addOption(Option.builder("pe").desc("Pair end mode (Optional)").build());
		options.addOption(
				Option.builder()
						.longOpt("format")
						.hasArg()
						.desc("specify format of input: mappedread, sam")
						.required()
						.build());

		CommandLine cmd = new CMDHelper(args, "cpmr [options]", options).build();

		String referenceChromosomeFileName = cmd.getOptionValue("r");
		String mappedReadFileName = cmd.getOptionValue("m");
		String outputPath = cmd.getOptionValue("o");
		double partial_methyl_threshold = Double.valueOf(cmd.getOptionValue("p"));
		boolean isPairEnd = cmd.hasOption("pe");
		int min_cpg_coverage = Integer.parseInt(cmd.getOptionValue("mcc"));
		int min_interval_cpg = Integer.parseInt(cmd.getOptionValue("mic"));
		int min_interval_reads = Integer.parseInt(cmd.getOptionValue("mir"));
		MappedReadFileFormat mappedReadFileFormat;
		switch (cmd.getOptionValue("format")) {
			case "mappedread":
				mappedReadFileFormat = MappedReadFileFormat.MAPPEDREAD;
				break;
			case "sam":
				mappedReadFileFormat = MappedReadFileFormat.SAM;
				break;
			default:
				throw new RuntimeException("unknown input format:\t" + cmd.getOptionValue("format"));
		}

		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = IOUtils.readReferenceChromosome(referenceChromosomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), MethylationUtils.REFERENCE_INIT_POS);
		System.out.printf("Load refChr %s-%d-%d with %d RefCpGs. Complete in %f s\n", refChr.getChr(),
				refChr.getStart(), refChr.getEnd(),
				refCpGList.size(), (System.currentTimeMillis() - start) / 1000.0);

		// load mapped reads
		start = System.currentTimeMillis();
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		List<MappedRead> mappedReadList = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
				new MappedReadLineProcessor(isPairEnd, mr -> mr.generateCpGsInRead(refMap) > 0, mappedReadFileFormat));
		System.out.printf("Load %d Mapped Reads. Complete in %f s\n", mappedReadList.size(),
				(System.currentTimeMillis() - start) / 1000.0);

		if (mappedReadList.size() == 0) {
			System.err.println(
					"0 Mapped Reads loaded! Please double check your reference and reads file. Make sure they use same reference chromosome name");
			System.exit(Constant.EXIT_ON_ERROR);
		}

		InputReadsSummary cpgReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		mappedReadList.forEach(cpgReadsSummary::addMappedRead);

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
		writeReport(reportFileName, cpgReadsSummary.getSummaryString(
				"\nSummary of reads with at least 1 CpG:\n") + intervalReadsSummary.getSummaryString(
				"\nSummary of reads in interval:\n") + InputReadsSummary.getCpGCoverageSummary(
				refCpGCollection),
				refCpGList.size(), cpmr.getRawIntervalCount(), immutableGenomicIntervals.size());
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
		// TODO mark CpG site on reference string, e.g. GaattcCGaattaC
		mappedReadWriter.write(String.format("ref:\t%s\n",
				immutableGenomicInterval.getRefString()
						.substring(immutableGenomicInterval.getStart() - MethylationUtils.REFERENCE_INIT_POS,
								immutableGenomicInterval.getEnd() + 1 - MethylationUtils.REFERENCE_INIT_POS)));
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
				// bed file is 0-based start and end.
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