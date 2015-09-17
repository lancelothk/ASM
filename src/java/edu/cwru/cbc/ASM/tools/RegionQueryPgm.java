package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.genomicInterval.MethylStatBedInterval;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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

		List<MethylStatBedInterval> regions = getMethylStatBedIntervals(bedFileName);

		List<String> lines = Files.readLines(new File(inputFileName), Charsets.UTF_8);
		lines.stream().forEach(line -> {
			List<String> itemList = tabSplitter.splitToList(line);
			String chr = itemList.get(0);
			int pos = Integer.parseInt(itemList.get(1));
			int coverage = Integer.parseInt(itemList.get(2));
			double methylLevel = Double.parseDouble(itemList.get(4));
			regions.forEach(r -> r.intersect(chr, pos, coverage, methylLevel));
		});
//				checkEachCpG(inputFileName, regions);

		Files.asCharSink(new File(outputFileName), Charsets.UTF_8).writeLines(regions.stream().map(
				MethylStatBedInterval::toString).collect(
				Collectors.toList()));
	}

	private static void checkEachCpG(String inputFileName, final List<MethylStatBedInterval> regions) throws
			IOException {
		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor() {

			@Override
			public boolean processLine(String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				String chr = itemList.get(0);
				int pos = Integer.parseInt(itemList.get(1));
				int coverage = Integer.parseInt(itemList.get(2));
				double methylLevel = Double.parseDouble(itemList.get(4));
				regions.forEach(r -> r.intersect(chr, pos, coverage, methylLevel));
				return true;
			}

			@Override
			public Object getResult() {
				return null;
			}
		});
	}

	private static List<MethylStatBedInterval> getMethylStatBedIntervals(String bedFileName) throws
			IOException {
		return Files.readLines(new File(bedFileName), Charsets.UTF_8,
				new LineProcessor<List<MethylStatBedInterval>>() {
					private List<MethylStatBedInterval> regions = new ArrayList<>();

					@Override
					public boolean processLine(String line) throws IOException {
						List<String> itemList = tabSplitter.splitToList(line);
						regions.add(new MethylStatBedInterval(itemList.get(0), Integer.parseInt(itemList.get(1)),
								Integer.parseInt(itemList.get(2)), itemList.get(3)));
						return true;
					}

					@Override
					public List<MethylStatBedInterval> getResult() {
						return regions;
					}
				});
	}
}