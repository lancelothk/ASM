package edu.cwru.cbc.ASM.commons.bed;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;

import javax.annotation.Nonnull;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by lancelothk on 4/2/15.
 * Utils for handling Bed format file.
 */
//TODO refactor & unit test
public class BedUtils {
	private static final Splitter tabSplitter = Splitter.on("\t");
	public static List<BedInterval> readSingleChromBedRegions(String bedFileName) throws IOException {
		Map<String, List<BedInterval>> regionsMap = readBedRegions(bedFileName);
		if (regionsMap.size() == 0) {
			return new ArrayList<>();
		} else if (regionsMap.size() == 1) {
			return Iterables.get(regionsMap.values(), 0);
		} else {
			throw new RuntimeException("Bed file contains regions from multiple chromosomes!");
		}
	}

	/**
	 * read Bed format regions.
	 * Only include first 4 columns: chr, start, end, name.
	 * And the last column label if readLabel is true.
	 */
	public static Map<String, List<BedInterval>> readBedRegions(String bedFileName) throws IOException {
		return Files.readLines(new File(bedFileName), Charsets.UTF_8,
				new LineProcessor<Map<String, List<BedInterval>>>() {
					private Map<String, List<BedInterval>> genomicIntervalMap = new HashMap<>();

					@Override
					public boolean processLine(@Nonnull String line) throws IOException {
						List<String> items = tabSplitter.splitToList(line);
						if (items.get(0).equals("chr") || line.equals("") || items.get(0).startsWith("#")) {
							// skip column name
							return true;
						}
						if (items.size() < 3) {
							throw new RuntimeException("less than 3 columns in bed file!");
						} else if (items.size() == 3) {
							// 3 columns bed format
							// suppose input bed is 0-based start and 1-based end.
							// internal representations of position are all 0-based.
							addRegionToList(new BedInterval(items.get(0), Integer.parseInt(items.get(1)),
									Integer.parseInt(items.get(2)) - 1, ""), genomicIntervalMap);
						} else if (items.size() == 4) {
							// 4 columns bed format
							// suppose input bed is 0-based start and 1-based end.
							// internal representations of position are all 0-based.
							addRegionToList(new BedInterval(items.get(0), Integer.parseInt(items.get(1)),
									Integer.parseInt(items.get(2)) - 1, items.get(3)), genomicIntervalMap);
						} else {
							// more than 4 columns bed format
							// suppose input bed is 0-based start and 1-based end.
							// internal representations of position are all 0-based.
							addRegionToList(new BedInterval(items.get(0), Integer.parseInt(items.get(1)),
									Integer.parseInt(items.get(2)) - 1, items.get(3),
									items.subList(4, items.size())), genomicIntervalMap);
						}
						return true;
					}

					@Override
					public Map<String, List<BedInterval>> getResult() {
						return genomicIntervalMap;
					}
				});
	}

	private static void addRegionToList(BedInterval newRegion,
	                                    Map<String, List<BedInterval>> genomicIntervalMap) {
		if (genomicIntervalMap.containsKey(newRegion.getChr())) {
			genomicIntervalMap.get(newRegion.getChr()).add(newRegion);
		} else {
			List<BedInterval> bedIntervalList = new ArrayList<>();
			bedIntervalList.add(newRegion);
			genomicIntervalMap.put(newRegion.getChr(), bedIntervalList);
		}
	}

	private static boolean hasOverlap(BedInterval regionA, BedInterval regionB) {
		return regionA.getChr().equals(regionB.getChr()) &&
				(regionA.getStart() <= regionB.getEnd() && regionA.getEnd() >= regionB.getStart());
	}

	public static Collection<BedInterval> intersect(Collection<BedInterval> regionsA,
	                                                Collection<BedInterval> regionsB) {
		List<BedInterval> intersections = new ArrayList<>();
		for (BedInterval regionA : regionsA) {
			for (BedInterval regionB : regionsB) {
				if (hasOverlap(regionA, regionB)) {
					regionA.addIntersectedRegion(regionB);
					regionB.addIntersectedRegion(regionA);
					int right = Math.min(regionA.getEnd(), regionB.getEnd());
					int left = Math.max(regionA.getStart(), regionB.getStart());
					intersections.add(new BedInterval(regionA.getChr(), left, right, regionA.getName() + "-" +
							regionB.getName()));
				}
			}
		}
		return intersections;
	}


	public static void writeBedRegions(List<BedInterval> regions, String bedFileName) throws IOException {
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(bedFileName));
		regions.sort(BedInterval::compareTo);
		for (BedInterval region : regions) {
			bedWriter.write(region.toBedString() + "\n");
		}
		bedWriter.close();
	}

	public static void writeBedWithIntersection(List<BedInterval> regionsList, String bedFileName) throws
			IOException {
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(bedFileName));
		for (BedInterval region : regionsList) {
			bedWriter.write(region.toBedWithIntersectionString() + "\n");
		}
		bedWriter.close();
	}
}
