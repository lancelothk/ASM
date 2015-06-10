package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.GenomicInterval.BedInterval;
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
public class RandomRegionSelection {
	public static void main(String[] args) throws IOException {
		String currUserHome = System.getProperty("user.home");
		String inputFile = currUserHome +
				"/experiments/ASM/simulation/CpGIslandsRegions/nonCpGIslandsRegions_selected_100000_15_5.Bed";
		int numberOfSelection = 100;

		Random rand = new Random();
		List<BedInterval> regionList = BedUtils.readSingleChromBedRegions(inputFile);
		List<BedInterval> selection = new ArrayList<>();
		while (selection.size() != numberOfSelection) {
			BedInterval bedInterval = regionList.get(rand.nextInt(regionList.size()));
			if (!selection.contains(bedInterval)) {
				selection.add(bedInterval);
			}
		}

		BufferedWriter writer = new BufferedWriter(
				new FileWriter(inputFile.replace(".Bed", "") + "_selection" + ".Bed"));
		for (BedInterval bedInterval : selection) {
			writer.write(bedInterval.toBedString() + "\n");
		}

		writer.close();
	}
}
