package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.Methylation.RefChr;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.util.BitSet;

/**
 * Created by lancelothk on 2/24/15.
 * Get the coverage details of a genomic region.
 */
public class RegionSelectionOnCoverage {
	private static int bpsPassedCoverageThreshold = 0;
	private static int totalLength = 0;
	private static int count = 0;
	private static int qualifiedCount = 0;

	public static void main(String[] args) throws IOException {
		String currUserHome = System.getProperty("user.home");
		RefChr refChr = CommonsUtils.readReferenceGenome(currUserHome + "/experiments/ASM/data/hg18_chr20.fa");
		String regionType = "nonCGI";
		String inputFolderName =
				currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/" + regionType + "_regions";
		String outputBedRegionFolderName = currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/";
		// >= min_coverage && < max_coverage
		execution(refChr, regionType, inputFolderName, outputBedRegionFolderName, 10, 5, 5);
	}

	private static void execution(RefChr refChr, String regionType, String inputFolderName, String outputBedRegionFolderName,
								  int max_coverage, int min_coverage, int min_cpg_number) throws IOException {
		bpsPassedCoverageThreshold = 0;
		totalLength = 0;
		count = 0;
		qualifiedCount = 0;

		System.out.printf("Coverage threshold:\t%d\t", min_coverage);

		String prefix = "i90_r1_chr20_";
		//		String inputFolderName = currUserHome + "/experiments/ASM/data/";
		//		String prefix = "i90_r1_chr20_";

		BufferedWriter writer = new BufferedWriter(new FileWriter(outputBedRegionFolderName +
				String.format(regionType + "Regions_selected_%d_%d_%d.Bed", max_coverage, min_coverage,
						min_cpg_number)));

		File inputFolder = new File(inputFolderName);
		if (!inputFolder.isDirectory()) {
			throw new RuntimeException("input is not a folder!");
		}

		File[] inputFiles = inputFolder.listFiles();
		assert inputFiles != null;
		for (File inputFile : inputFiles) {
			if (inputFile.getName().startsWith(prefix) && !inputFile.getName().endsWith(".aligned")) {
				//System.out.println(inputFile.getName());
				selectHighCoverageRegion(refChr, inputFile, writer, prefix, max_coverage, min_coverage, min_cpg_number);
			}
		}

		writer.close();
		System.out.printf("total count:%d\tqualifiedCount:%d\n", count, qualifiedCount);
		System.out.printf("bps passed coverage threshold:%d\ttotal length:%d\t%.2f\n", bpsPassedCoverageThreshold,
				totalLength, bpsPassedCoverageThreshold / (double) totalLength);
	}

	private static void selectHighCoverageRegion(RefChr refChr, File inputFile, BufferedWriter writer, String prefix,
			int max_coverage, int min_coverage, int min_cpg_number) throws IOException {
		String[] items = inputFile.getName().replace(prefix, "").split("-");
		if (items.length != 2) {
			throw new RuntimeException("invalid file name format!");
		}

		int regionStart = Integer.valueOf(items[0]);
		int regionEnd = Integer.valueOf(items[1]);


		int[] coverage = new int[regionEnd - regionStart + 1];

		totalLength += coverage.length;

		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		String line;
		while ((line = reader.readLine()) != null) {
			try {
				items = line.split("\t");
				int start = Integer.valueOf(items[2]);
				int end = Integer.valueOf(items[3]);
				for (int i = start; i <= end; i++) {
					if (i >= regionStart && i <= regionEnd) {
						coverage[i - regionStart]++;
					}
				}
			} catch (Exception e) {
				System.out.println(inputFile.getName());
				throw e;
			}
		}

		BitSet bitSet = new BitSet(regionEnd - regionStart + 1);
		for (int i = 0; i < coverage.length; i++) {
			if (coverage[i] >= min_coverage && coverage[i] < max_coverage) {
				bitSet.set(i);
				bpsPassedCoverageThreshold++;
			}
		}

		int nextTrue;
		int nextFalse = -1;
		while ((nextTrue = bitSet.nextSetBit(nextFalse + 1)) != -1) {
			if ((nextFalse = bitSet.nextClearBit(nextTrue + 1)) != -1) {
				int start = nextTrue + regionStart;
				int end = nextFalse - 1 + regionStart;
				int cpgCount = StringUtils.countMatches(refChr.getRefString().substring(start, end), "CG");
				if (cpgCount >= min_cpg_number) {
					System.out.printf("region:%d-%d\n", start, end);
					writer.write(String.format("%s\t%d\t%d\t%s\n", refChr.getChr(), start, end, cpgCount));
					qualifiedCount++;
				}
				count++;
			} else {
				int cpgCount = StringUtils.countMatches(refChr.getRefString()
								.substring(regionStart + nextTrue + regionStart,
										regionStart + nextFalse - 1 + regionStart), "CG");
				if (cpgCount >= min_cpg_number) {
					System.out.printf("region:%d-%d\n", nextTrue, regionEnd);
					writer.write(String.format("%s\t%d\t%d\t%s\n", refChr.getChr(), nextTrue, regionEnd, cpgCount));
					qualifiedCount++;
				}
				count++;
				break;
			}
		}
		reader.close();
	}
}
