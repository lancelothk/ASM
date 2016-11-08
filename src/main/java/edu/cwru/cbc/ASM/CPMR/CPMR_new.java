package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.genomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadHandler;
import edu.cwru.cbc.ASM.commons.io.MappedReadReader;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * main entry of CPMR program.
 */
public class CPMR_new {
	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("r").hasArg().desc("Reference Path").required().build());
		options.addOption(Option.builder("i").hasArg().desc("bam File").required().build());
		options.addOption(Option.builder("c").hasArg().desc("chromosome").required().build());
		options.addOption(Option.builder("o").hasArg().desc("Output Path").required().build());
		options.addOption(Option.builder("p").hasArg().desc("Partial methylation threshold").required().build());
		options.addOption(Option.builder("mcc").hasArg().desc("Minimum adjacent CpG coverage").required().build());
		options.addOption(Option.builder("mic").hasArg().desc("Minimum interval CpG number").required().build());
		options.addOption(Option.builder("mir").hasArg().desc("Minimum interval read number").required().build());
		options.addOption(Option.builder("pe").desc("pair end mode").build());
		options.addOption(
				Option.builder()
						.longOpt("format")
						.hasArg()
						.desc("specify format of input: mappedread, sam")
						.required()
						.build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String referencePath = cmd.getOptionValue("r");
		String chr = cmd.getOptionValue("c");
		String mappedReadFileName = cmd.getOptionValue("i");
		String outputPath = cmd.getOptionValue("o");
		double partial_methyl_threshold = Double.valueOf(cmd.getOptionValue("p"));
		if (partial_methyl_threshold < 0 || partial_methyl_threshold > 1) {
			System.err.println("invalid partial methylation threshold. Should be within [0,1]");
			System.exit(1);
		}
		boolean isPairEnd = cmd.hasOption("pe");
		int min_cpg_coverage = Integer.parseInt(cmd.getOptionValue("mcc"));
		if (min_cpg_coverage < 0) {
			System.err.println("invalid Minimum adjacent CpG coverage. Should be larger than zero");
			System.exit(1);
		}
		int min_interval_cpg = Integer.parseInt(cmd.getOptionValue("mic"));
		if (min_interval_cpg < 0) {
			System.err.println("invalid Minimum interval CpG number. Should be larger than zero");
			System.exit(1);
		}
		int min_interval_reads = Integer.parseInt(cmd.getOptionValue("mir"));
		if (min_interval_reads < 0) {
			System.err.println("invalid Minimum interval read number. Should be larger than zero");
			System.exit(1);
		}

		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = IOUtils.readReferenceChromosome(referencePath + "/" + chr + ".fa");
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), MethylationUtils.REFERENCE_INIT_POS);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		// load mapped reads
		start = System.currentTimeMillis();
		String summaryFileName = outputPath + "/CPMR.bed";
		String reportFileName = outputPath + "/CPMR.report";
		File intervalFolder = new File(outputPath);
		if (!intervalFolder.exists()) {
			if (!intervalFolder.mkdirs()) {
				throw new RuntimeException("failed to create folder:\t" + intervalFolder);
			}
		}
		if (refCpGList.size() < 2) {
			System.err.println("reference contains less than three CpG sites!");
			System.exit(1);
		}
		Counter counter = new Counter();
		InputReadsSummary cpgReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		InputReadsSummary intervalReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		List<RefCpG> refCpGCollection = new ArrayList<>();
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		intervalSummaryWriter.write("chr\tstart\tend\treadCount\tCpGCount\n");
		// require input file is sorted.
		MappedReadReader.readMappedReads(new File(mappedReadFileName), refMap, new MappedReadHandler() {
			int refCpGIndex = 0;
			List<MappedRead> mappedReads = new ArrayList<>();
			List<RefCpG> intervalRefCpGs = new ArrayList<>();

			@Override
			public void processMappedRead(MappedRead mappedRead) throws IOException {
				cpgReadsSummary.addMappedRead(mappedRead);
				// stop at last RefCpG
				if (refCpGIndex == refCpGList.size() - 1) {
					writeQualifiedInterval();
					return;
				}
				// if read doesn't cover curr or next RefCpG
				while (mappedRead.getStart() > refCpGList.get(refCpGIndex).getPos()
						&& mappedRead.getStart() > refCpGList.get(refCpGIndex + 1).getPos()) {
					refCpGIndex++;
				}
				// read cover curr and next AND Consecutively Partially methylated(CPM)
				if (mappedRead.getStart() < refCpGList.get(refCpGIndex).getPos()
						&& mappedRead.getStart() < refCpGList.get(refCpGIndex + 1).getPos()
						&& refCpGList.get(refCpGIndex).hasCommonRead(refCpGList.get(refCpGIndex + 1), min_cpg_coverage)
						&& refCpGList.get(refCpGIndex).hasPartialMethyl(partial_methyl_threshold)
						&& refCpGList.get(refCpGIndex + 1).hasPartialMethyl(partial_methyl_threshold)) {
					mappedReads.add(mappedRead);
				} else {
					// cover (curr or next) OR not qualified Consecutively Partially methylated
					writeQualifiedInterval();
					refCpGIndex++;
					mappedReads = new ArrayList<>();
					mappedReads.add(mappedRead);
					intervalRefCpGs = new ArrayList<>();
					intervalRefCpGs.add(refCpGList.get(refCpGIndex + 1));
				}
			}

			private void writeQualifiedInterval() throws IOException {
				counter.rawIntervalCount++;
				if (intervalRefCpGs.size() >= min_interval_cpg && mappedReads.size() >= min_interval_reads) {
					counter.intervalCount++;
					// bed file is 0-based start and end.
					ImmutableGenomicInterval interval = new ImmutableGenomicInterval(refChr.getChr(),
							refChr.getRefString(),
							intervalRefCpGs.get(0).getPos(),
							intervalRefCpGs.get(intervalRefCpGs.size() - 1).getPos(), intervalRefCpGs,
							mappedReads);
					intervalSummaryWriter.write(
							String.format("%s\t%d\t%d\t%d\t%d\n", interval.getChr(), interval.getStart(),
									interval.getEnd(), interval.getMappedReadList().size(),
									interval.getRefCpGList().size()));
					writeMappedReadInInterval(outputPath, interval);
					mappedReads.forEach(intervalReadsSummary::addMappedRead);
					refCpGCollection.addAll(intervalRefCpGs);
				}
			}
		}, chr);
		intervalSummaryWriter.close();

		writeReport(reportFileName, cpgReadsSummary.getSummaryString(
				"\nSummary of reads with at least 1 CpG:\n") + intervalReadsSummary.getSummaryString(
				"\nSummary of reads in interval:\n") + InputReadsSummary.getCpGCoverageSummary(
				refCpGCollection),
				refCpGList.size(), counter.rawIntervalCount, counter.intervalCount);
		System.out.println("CPMR complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");
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

	static class Counter {
		int rawIntervalCount;
		int intervalCount;
	}

}