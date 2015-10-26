package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.genomicInterval.MethylStatBedInterval;
import edu.cwru.cbc.ASM.commons.methylation.RefCpGStat;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by kehu on 9/16/15.
 * Similar to bigWigAverageOverBed, get stat over a bed region.
 */
public class RegionQueryPgm {
	private static final Splitter tabSplitter = Splitter.on("\t");

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input File").required().build());
		options.addOption(Option.builder("b").hasArg().desc("bed file").required().build());
		options.addOption(Option.builder("o").hasArg().desc("output File").required().build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String inputFileName = cmd.getOptionValue("i");
		String bedFileName = cmd.getOptionValue("b");
		String outputFileName = cmd.getOptionValue("o");

		Map<String, List<MethylStatBedInterval>> regionMap = getMethylStatBedIntervals(bedFileName);

		Map<String, HashIntObjMap<RefCpGStat>> refMap = getRefCpGStatHashIntObjMap(inputFileName);

		Files.asCharSink(new File(outputFileName), Charsets.UTF_8).writeLines(
				regionMap.entrySet()
						.stream()
						.filter(entry -> refMap.containsKey(entry.getKey()))
						.peek(entry -> entry.getValue()
								.stream()
								.forEach(r -> r.queryRefMap(refMap.get(entry.getKey()))))
						.map(Map.Entry::getValue)
						.flatMap(Collection::stream)
						.map(MethylStatBedInterval::toString)
						.collect(Collectors.toList()));
	}

	private static Map<String, HashIntObjMap<RefCpGStat>> getRefCpGStatHashIntObjMap(String inputFileName) throws
			IOException {
		return Files.readLines(new File(inputFileName), Charsets.UTF_8,
				new LineProcessor<Map<String, HashIntObjMap<RefCpGStat>>>() {
					private Map<String, HashIntObjMap<RefCpGStat>> refMap = new HashMap<>();

					@Override
					public boolean processLine(@Nonnull String line) throws IOException {
						List<String> itemList = tabSplitter.splitToList(line);
						String chr = itemList.get(0);
						int pos = Integer.parseInt(itemList.get(1));
						int coverage = Integer.parseInt(itemList.get(2));
						double methylLevel = Double.parseDouble(itemList.get(4));
						HashIntObjMap<RefCpGStat> refCpGMap = refMap.get(chr);
						if (refCpGMap == null) {
							refCpGMap = HashIntObjMaps.newMutableMap();
							refCpGMap.put(pos, new RefCpGStat(pos, coverage, methylLevel));
							refMap.put(chr, refCpGMap);
						} else {
							refCpGMap.put(pos, new RefCpGStat(pos, coverage, methylLevel));
						}
						return true;
					}

					@Override
					public Map<String, HashIntObjMap<RefCpGStat>> getResult() {
						return refMap;
					}
				});
	}

	private static Map<String, List<MethylStatBedInterval>> getMethylStatBedIntervals(String bedFileName) throws
			IOException {
		return Files.readLines(new File(bedFileName), Charsets.UTF_8,
				new LineProcessor<Map<String, List<MethylStatBedInterval>>>() {
					private Map<String, List<MethylStatBedInterval>> regionMap = new HashMap<>();

					@Override
					public boolean processLine(@Nonnull String line) throws IOException {
						List<String> itemList = tabSplitter.splitToList(line);
						String chr = itemList.get(0);
						List<MethylStatBedInterval> regionList = regionMap.get(chr);
						if (regionList == null) {
							regionList = new ArrayList<>();
							regionList.add(new MethylStatBedInterval(itemList.get(0), Integer.parseInt(itemList.get(1)),
									Integer.parseInt(itemList.get(2)), itemList.get(3)));
							regionMap.put(chr, regionList);
						} else {
							regionList.add(new MethylStatBedInterval(itemList.get(0), Integer.parseInt(itemList.get(1)),
									Integer.parseInt(itemList.get(2)), itemList.get(3)));
						}
						return true;
					}

					@Override
					public Map<String, List<MethylStatBedInterval>> getResult() {
						return regionMap;
					}
				});
	}
}