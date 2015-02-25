package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.DataType.GenomicRegion;
import edu.cwru.cbc.ASM.commons.Utils;

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
		String inputFile = "/home/lancelothk/experiments/ASM/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength0.6.bed";

		Random rand = new Random();
		List<GenomicRegion> regionList = Utils.readBedRegions(inputFile);
		int numberOfSelection = 100;
		List<GenomicRegion> selection = new ArrayList<>();
		for (int i = 0; i < numberOfSelection; i++) {
			selection.add(regionList.get(rand.nextInt(regionList.size())));
		}


		BufferedWriter writer = new BufferedWriter(new FileWriter(inputFile.replace(".bed", "") + "_selection" + ".bed"));
		for (GenomicRegion genomicRegion : selection) {
			writer.write(genomicRegion.toBedString() + "\n");
		}

		writer.close();
	}
}
