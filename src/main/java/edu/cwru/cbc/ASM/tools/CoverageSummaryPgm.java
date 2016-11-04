package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.CPMR.InputReadsSummary;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadReader;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by kehu on 11/4/16.
 * Calculate coverage summary of given SAM/BAM file.
 */
public class CoverageSummaryPgm {
	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("r").hasArg().desc("Reference File Path").required().build());
		options.addOption(Option.builder("i").hasArg().desc("Input SAM/BAM file").required().build());
		options.addOption(
				Option.builder("p").hasArg().desc("Genomic position, e.g. chr1, chr1:1000-2000.").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceFilePath = cmd.getOptionValue("r");
		String inputFilePath = cmd.getOptionValue("i");

		String chr;
		int start = -1, end = -1;
		boolean hasPos = false;
		String[] items = cmd.getOptionValue("p").split(":|-");
		if (items.length == 0) {
			System.err.println("invalid position!");
			System.exit(1);
		}
		chr = items[0];
		if (items.length == 3) {
			start = Integer.parseInt(items[1]);
			end = Integer.parseInt(items[2]);
			if ((start < 0 && end >= 0) || (start >= 0 && end < 0)) {
				System.err.println("invalid position!");
				System.exit(1);
			} else {
				hasPos = true;
			}
		} else if (items.length != 1) {
			System.err.println("invalid position!");
			System.exit(1);

		}

		long startTime = System.currentTimeMillis();
		RefChr refChr = IOUtils.readReferenceChromosome(referenceFilePath);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), MethylationUtils.REFERENCE_INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - startTime) / 1000.0 + "s");


		startTime = System.currentTimeMillis();
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		List<MappedRead> mappedReads;
		if (!hasPos) {
			mappedReads = MappedReadReader.readMappedReads(new File(inputFilePath), chr);
		} else {
			mappedReads = MappedReadReader.readMappedReads(new File(inputFilePath), chr, start, end);
		}

		System.out.println(
				"load mappedReadLinkedHashMap complete\t" + (System.currentTimeMillis() - startTime) / 1000.0 + "s");

		InputReadsSummary cpgReadsSummary = new InputReadsSummary(refChr.getRefString().length());
		mappedReads.forEach(cpgReadsSummary::addMappedRead);
		System.out.println(cpgReadsSummary.getSummaryString());
	}
}
