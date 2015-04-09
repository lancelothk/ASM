package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.GenomicInterval;
import edu.cwru.cbc.ASM.commons.bed.BedUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by lancelothk on 2/24/15.
 * random select regions.
 */
public class RandomRegions {
	public static void main(String[] args) throws IOException {
		String currUserHome = System.getProperty("user.home");
		String inputFile = currUserHome +
				"/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength_8_0.2.Bed";

		Random rand = new Random();
		List<GenomicInterval> regionList = BedUtils.readSingleChromBedRegions(inputFile);
		int numberOfSelection = 100;
		List<GenomicInterval> selection = new ArrayList<>();
		while (selection.size() != numberOfSelection) {
			GenomicInterval genomicInterval = regionList.get(rand.nextInt(regionList.size()));
			if (!selection.contains(genomicInterval)) {
				selection.add(genomicInterval);
			}
		}


		BufferedWriter writer = new BufferedWriter(
				new FileWriter(inputFile.replace(".Bed", "") + "_selection" + ".Bed"));
		for (GenomicInterval genomicInterval : selection) {
			writer.write(genomicInterval.toBedString() + "\n");
		}

		writer.close();
	}
}
