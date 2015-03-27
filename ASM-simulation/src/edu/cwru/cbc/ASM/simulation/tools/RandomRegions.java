package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.DataType.GenomicInterval;

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
                "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength_8_0.2.bed";

		Random rand = new Random();
		List<GenomicInterval> regionList = CommonsUtils.readBedRegions(inputFile);
		int numberOfSelection = 100;
		List<GenomicInterval> selection = new ArrayList<>();
		for (int i = 0; selection.size() != numberOfSelection; i++) {
			GenomicInterval genomicInterval = regionList.get(rand.nextInt(regionList.size()));
			if (!selection.contains(genomicInterval)) {
				selection.add(genomicInterval);
			}
		}


		BufferedWriter writer = new BufferedWriter(new FileWriter(inputFile.replace(".bed", "") + "_selection" + ".bed"));
		for (GenomicInterval genomicInterval : selection) {
			writer.write(genomicInterval.toBedString() + "\n");
		}

		writer.close();
	}
}
