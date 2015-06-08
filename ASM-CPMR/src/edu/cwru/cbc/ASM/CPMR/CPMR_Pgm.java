package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.Read.MappedReadLineProcessorWithSummary;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * main entry of program.
 */
public class CPMR_Pgm {
	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption("r", true, "Reference File");
		options.addOption("m", true, "MappedRead File");
		options.addOption("o", true, "Output Path");
		options.addOption("mcc", true, "Minimum adjacent CpG coverage");
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


		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = CommonsUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), CPMR.INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		// load mapped reads
		start = System.currentTimeMillis();
		String reportString = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
				new MappedReadLineProcessorWithSummary(refCpGList, min_read_cpg, refChr.getRefString().length()));
		System.out.println("load mappedReadList complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		String summaryFileName = outputPath + "/CPMR.bed";
		String reportFileName = outputPath + "/CPMR.report";
		File intervalFolder = new File(outputPath);
		if (!intervalFolder.exists()) {
			//noinspection ResultOfMethodCallIgnored
			intervalFolder.mkdirs();
		}

		CPMR cpmr = new CPMR(refChr, refCpGList, min_cpg_coverage, min_read_cpg,
				min_interval_cpg, min_interval_reads);
		List<List<RefCpG>> cpgSiteIntervalList = cpmr.getIntervals();
		cpmr.writeIntervals(outputPath, summaryFileName, cpgSiteIntervalList);
		cpmr.writeReport(reportFileName, reportString, cpgSiteIntervalList);
	}

}