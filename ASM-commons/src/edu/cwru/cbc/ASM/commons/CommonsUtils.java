package edu.cwru.cbc.ASM.commons;

import com.google.common.base.Charsets;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.GenomicInterval;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Created by kehu on 11/13/14.
 * CommonsUtils for methods shared by modules under ASM project.
 */
public class CommonsUtils {

	/**
	 * Extract reference CpG from reference string.
	 *
	 * @param reference should be upperCase string and without space in it.
	 * @return returned List is sorted by position.
	 */
	public static List<RefCpG> extractCpGSite(String reference, int initPos) {
		// initialize refCpGList with 1/20 of reference size, since the probability of CG occurrence should be 1/16.
		// Actual probability is higher than 1/16 from current observation. In CPG dense region, it should be much higher.
		List<RefCpG> reFCpGList = new ArrayList<>(reference.length() / 20);
		for (int i = 0; i < reference.length() - 1; i++) {
			if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
				reFCpGList.add(new RefCpG(i + initPos));
			}
		}
		return reFCpGList;
	}

	/**
	 * Read Reference genome into RefChr type.
	 * The ref string is concerted to upper case.
	 *
	 * @throws IOException
	 */
	public static RefChr readReferenceGenome(String inputFileName) throws IOException {
		return Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor<RefChr>() {
			private String chr;
			private StringBuilder referenceBuilder = new StringBuilder();

			@Override
			public boolean processLine(String s) throws IOException {
				if (s.startsWith(">")) {
					chr = s.substring(1, s.length());
				} else {
					if (Pattern.compile("[^ACGTNacgtn]").matcher(s).find()) {
						throw new RuntimeException("invalid character in sequence! only acgtn/ACGTN are allowed!");
					}
					referenceBuilder.append(s.replaceAll(" ", "").toUpperCase());
				}
				return true;
			}

			@Override
			public RefChr getResult() {
				if (chr == null) {
					throw new RuntimeException("didn't find reference name!");
				}
				return new RefChr(chr, referenceBuilder.toString());
			}
		});
	}

	public static List<GenomicInterval> readSingleChromBedRegions(String bedFileName) throws IOException {
		return readSingleChromBedRegions(bedFileName, false);
	}

	public static List<GenomicInterval> readSingleChromBedRegions(String bedFileName,
																  boolean readLabel) throws IOException {
		Map<String, List<GenomicInterval>> regionsMap = readBedRegions(bedFileName, readLabel);
		if (regionsMap.size() != 1) {
			throw new RuntimeException("bed file should only contain regions from single chromosome!");
		}
		return Iterables.get(regionsMap.values(), 0);
	}

	/**
	 * read bed format regions.
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

	/**
	 * read bed format regions.
	 * Only include first 4 columns: chr, start, end, name.
	 */
	public static Map<String, List<GenomicInterval>> readBedRegions(String bedFileName) throws IOException {
		return readBedRegions(bedFileName, false);
	}
}
