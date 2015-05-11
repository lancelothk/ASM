package edu.cwru.cbc.ASM.commons.bed;

import com.google.common.base.Charsets;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.GenomicInterval;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by lancelothk on 4/2/15.
 * Utils for handling Bed format file.
 */
public class BedUtils {
	public static List<GenomicInterval> readSingleChromBedRegions(String bedFileName) throws IOException {
		return readSingleChromBedRegions(bedFileName, false);
	}

	public static List<GenomicInterval> readSingleChromBedRegions(String bedFileName,
			boolean readLabel) throws IOException {
		Map<String, List<GenomicInterval>> regionsMap = readBedRegions(bedFileName, readLabel);
		if (regionsMap.size() != 1) {
			throw new RuntimeException("Bed file should only contain regions from single chromosome!");
		}
		return Iterables.get(regionsMap.values(), 0);
	}

	/**
	 * read Bed format regions.
	 * Only include first 4 columns: chr, start, end, name.
	 * And the last column label if readLabel is true.
	 */
	public static Map<String, List<GenomicInterval>> readBedRegions(String bedFileName,
			boolean readLabel) throws IOException {
		return Files.readLines(new File(bedFileName), Charsets.UTF_8,
				new LineProcessor<Map<String, List<GenomicInterval>>>() {
					private Map<String, List<GenomicInterval>> genomicIntervalMap = new HashMap<>();

					@Override
					public boolean processLine(String line) throws IOException {
						String[] items = line.split("\t");
						if (items[0].equals("chr") || line.equals("")) {
							// skip column name
							return true;
						}
						if (items.length > 4 && readLabel) {
							boolean isPositive;
							switch (line.charAt(line.length() - 1)) {
								case '+':
									isPositive = true;
									break;
								case '-':
									isPositive = false;
									break;
								case '*':
									isPositive = false;
									break;
								default:
									throw new RuntimeException("invalid label!\t" + line);
							}
							addRegionToList(new GenomicInterval(items[0].replace("chr", ""), Integer.parseInt(items[1]),
									Integer.parseInt(items[2]), items[3], isPositive), genomicIntervalMap);
						} else {
							addRegionToList(new GenomicInterval(items[0].replace("chr", ""), Integer.parseInt(items[1]),
									Integer.parseInt(items[2]), items[3]), genomicIntervalMap);
						}
						return true;
					}

					@Override
					public Map<String, List<GenomicInterval>> getResult() {
						return genomicIntervalMap;
					}
				});
	}

	private static void addRegionToList(GenomicInterval newRegion,
			Map<String, List<GenomicInterval>> genomicIntervalMap) {
		if (genomicIntervalMap.containsKey(newRegion.getChr())) {
			if (!hasOverlap(newRegion, genomicIntervalMap)) {
				genomicIntervalMap.get(newRegion.getChr()).add(newRegion);
			} else {
				throw new RuntimeException("Regions have overlap!");
			}
		} else {
			List<GenomicInterval> genomicIntervalList = new ArrayList<>();
			genomicIntervalList.add(newRegion);
			genomicIntervalMap.put(newRegion.getChr(), genomicIntervalList);
		}
	}

	private static boolean hasOverlap(GenomicInterval region, Map<String, List<GenomicInterval>> genomicIntervalMap) {
		return genomicIntervalMap.get(region.getChr()).stream().anyMatch(r -> hasOverlap(region, r));
	}

	/**
	 * Check if there exists any overlap between given region and regions
	 */
	private static boolean hasOverlap(GenomicInterval region, Collection<GenomicInterval> regions) {
		return regions.stream().anyMatch(r -> hasOverlap(region, r));
	}

	/**
	 * Check if there exists overlap in given regions.
	 */
	private static boolean hasOverlap(Collection<GenomicInterval> regions) {
		return regions.stream().anyMatch(r -> hasOverlap(r, regions));
	}

	private static boolean hasOverlap(GenomicInterval regionA, GenomicInterval regionB) {
		return regionA.getChr().equals(regionB.getChr()) &&
				(regionA.getStart() <= regionB.getEnd() && regionA.getEnd() >= regionB.getStart());
	}

	public static Collection<GenomicInterval> intersect(Collection<GenomicInterval> regionsA,
			Collection<GenomicInterval> regionsB) {
		List<GenomicInterval> intersections = new ArrayList<>();
		for (GenomicInterval regionA : regionsA) {
			for (GenomicInterval regionB : regionsB) {
				if (hasOverlap(regionA, regionB)) {
					regionA.addIntersectedRegion(regionB);
					regionB.addIntersectedRegion(regionA);
					int right = Math.min(regionA.getEnd(), regionB.getEnd());
					int left = Math.max(regionA.getStart(), regionB.getStart());
					intersections.add(new GenomicInterval(regionA.getChr(), left, right, regionA.getName() + "-" +
							regionB.getName()));
				}
			}
		}
		return intersections;
	}


	public static void writeBedRegions(Collection<GenomicInterval> regions, String bedFileName) throws IOException {
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(bedFileName));
		for (GenomicInterval region : regions) {
			bedWriter.write(region.toBedString() + "\n");
		}
		bedWriter.close();
	}

	public static void writeBedWithIntersection(Collection<GenomicInterval> regions, String bedFileName) throws IOException {
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(bedFileName));
		for (GenomicInterval region : regions) {
			bedWriter.write(region.toBedWithIntersectionString() + "\n");
		}
		bedWriter.close();
	}
}
