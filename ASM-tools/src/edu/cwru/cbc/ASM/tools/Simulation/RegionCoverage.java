package edu.cwru.cbc.ASM.tools.Simulation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by lancelothk on 2/24/15.
 * Get the coverage details of a genomic region.
 */
public class RegionCoverage {
	public static void main(String[] args) throws IOException {
		String currUserHome = System.getProperty("user.home");
		String inputFolderName = currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/";
		String prefix = "i90_r1_chr20_";

		File inputFolder = new File(inputFolderName);
		if (!inputFolder.isDirectory()){
			throw new RuntimeException("input is not a folder!");
		}

		File[] inputFiles = inputFolder.listFiles();
		for (File inputFile : inputFiles) {
			if (inputFile.getName().startsWith(prefix) && !inputFile.getName().endsWith(".aligned")){
				printQualifiedLength(inputFile, prefix);
			}
		}


	}

	private static void printQualifiedLength(File inputFile, String prefix) throws IOException {
		String[] items = inputFile.getName().replace(prefix, "").split("-");
		if (items.length != 2){
			throw new RuntimeException("invalid file name format!");
		}

		int regionStart = Integer.valueOf(items[0]);
		int regionEnd = Integer.valueOf(items[1]);


		int[] coverage = new int[regionEnd - regionStart + 1];

		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		String line;
		while ((line = reader.readLine()) != null) {
			items = line.split("\t");
			int start = Integer.valueOf(items[2]);
			int end = Integer.valueOf(items[3]);
			for (int i = start; i <= end; i++) {
				if (i >= regionStart && i <= regionEnd) {
					coverage[i - regionStart]++;
				}
			}
		}

		int count = 0;
		for (int i : coverage) {
			//System.out.print(i + "\t");
            if (i >= 8) {
                count++;
			}
		}
		//System.out.println();
		System.out.println(inputFile.getName() + "\t" + count);

		reader.close();
	}
}
