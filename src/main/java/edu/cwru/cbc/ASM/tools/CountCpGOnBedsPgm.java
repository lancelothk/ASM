package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by kehu on 7/27/16.
 * Count number of CpG sites on bed regions
 */
public class CountCpGOnBedsPgm {
	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption(Option.builder("r").hasArg().desc("Reference File Path").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceFilePath = cmd.getOptionValue("r");
		long startTime = System.currentTimeMillis();
		Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap = MethylationUtils.initializeGenomeRefCpGMap(
				IOUtils.readReferenceGenome(referenceFilePath));
		System.out.println("read reference finished: " + (System.currentTimeMillis() - startTime) / 1000);
		for (String bedFileName : cmd.getArgList()) {
			System.out.println(bedFileName);
			bedFileCpGCount(genomeRefCpGMap, bedFileName);
		}
	}

	private static void bedFileCpGCount(Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap, String bedFileName) throws
			IOException {
		String intervalCpGCountFile = bedFileName + ".CpGCount";
		Map<String, List<BedInterval>> bedRegionMap = BedUtils.readBedRegions(bedFileName);
		for (List<BedInterval> bedIntervals : bedRegionMap.values()) {
			for (BedInterval bedInterval : bedIntervals) {
				for (int i = bedInterval.getStart(); i <= bedInterval.getEnd(); i++) {
					RefCpG refCpG = genomeRefCpGMap.get(bedInterval.getChr()).get(i);
					if (refCpG != null) {
						bedInterval.addRefCpG(refCpG);
					}
				}
			}
		}
		List<BedInterval> allIntervalList = bedRegionMap.values()
				.stream()
				.flatMap(Collection::stream)
				.sorted()
				.collect(Collectors.toList());
		BufferedWriter writer = new BufferedWriter(new FileWriter(intervalCpGCountFile));
		for (BedInterval bedInterval : allIntervalList) {
			writer.write(String.format("%s\t%d\t%d\t%d\n", bedInterval.getChr(), bedInterval.getStart(),
					bedInterval.getEnd() + 1,
					bedInterval.getRefCpGList().size()));// internal bed is 0-based. Output bed end is 1-based.
		}
		System.out.println("median length: " + allIntervalList.stream()
				.mapToInt(i -> i.getEnd() - i.getStart() + 1).sorted()
				.toArray()[allIntervalList.size() / 2]);
		System.out.println("median CpG number: " + allIntervalList.stream()
				.mapToInt(i -> i.getRefCpGList().size()).sorted()
				.toArray()[allIntervalList.size() / 2]);
		writer.close();
	}
}
