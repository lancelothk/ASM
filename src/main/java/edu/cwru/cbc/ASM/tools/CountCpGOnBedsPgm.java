package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;

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
	public static void main(String[] args) throws IOException {
		String referenceFilePath = "/home/kehu/dataset/hg19";
		String bedFileName = "/home/kehu/experiments/ASM/result/ASM/whole/union.bed.whole";
		String intervalCpGCountFile = "/home/kehu/bedRegionWithCpGCount";
		Map<String, List<BedInterval>> bedRegionMap = BedUtils.readBedRegions(bedFileName);
		Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap = MethylationUtils.initializeGenomeRefCpGMap(
				IOUtils.readReferenceGenome(referenceFilePath));
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
					bedInterval.getEnd(), bedInterval.getRefCpGList().size()));
		}
		writer.close();
	}
}
