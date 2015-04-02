package edu.cwru.cbc.ASM.commons.Bed;

import com.google.common.base.Charsets;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.GenomicInterval;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
						if (readLabel) {
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
							addRegionToList(new GenomicInterval(items[0], Integer.parseInt(items[1]),
									Integer.parseInt(items[2]), items[3], isPositive), genomicIntervalMap);
						} else {
							addRegionToList(new GenomicInterval(items[0], Integer.parseInt(items[1]),
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
		return genomicIntervalMap.get(region.getChr())
				.stream()
				.anyMatch(r -> region.getStart() <= r.getEnd() && region.getEnd() >= r.getStart());
	}

	//	/**
	//	 * Check if there exists any overlap between given region and regions
	//	 */
	//	private static boolean hasOverlap(GenomicInterval region, Collection<GenomicInterval> regions) {
	//		return regions.stream()
	//				.anyMatch(r -> region.getStart() <= r.getEnd() && region.getEnd() >= r.getStart());
	//	}
	//
	//	/**
	//	 * Check if there exists overlap in given regions.
	//	 */
	//	private static boolean hasOverlap(Collection<GenomicInterval> regions) {
	//		return regions.stream().anyMatch(r-> hasOverlap(r, regions));
	//	}

}
