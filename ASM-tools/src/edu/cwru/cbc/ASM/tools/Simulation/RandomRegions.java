package edu.cwru.cbc.ASM.tools.Simulation;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.DataType.GenomicRegion;

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
		List<GenomicRegion> regionList = CommonsUtils.readBedRegions(inputFile);
		int numberOfSelection = 100;
		List<GenomicRegion> selection = new ArrayList<>();
		for (int i = 0; selection.size() != numberOfSelection; i++) {
			GenomicRegion genomicRegion = regionList.get(rand.nextInt(regionList.size()));
			if (!selection.contains(genomicRegion)) {
				selection.add(genomicRegion);
			}
		}


		BufferedWriter writer = new BufferedWriter(new FileWriter(inputFile.replace(".bed", "") + "_selection" + ".bed"));
		for (GenomicRegion genomicRegion : selection) {
			writer.write(genomicRegion.toBedString() + "\n");
		}

		writer.close();
	}
}
