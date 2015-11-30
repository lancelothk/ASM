package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by kehu on 7/6/15.
 * Tool to merge adjacent bed regions.
 */
public class MergeBedRegion {

	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input bed file").build());
		options.addOption(Option.builder("r").hasArg().desc("range to merge").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String bedFileName = cmd.getOptionValue("i");
		int range = Integer.parseInt(cmd.getOptionValue("r"));
		execute(range, bedFileName);
	}

	private static void execute(int range, String bedFileName) throws IOException {
		Map<String, List<BedInterval>> bedRegionMap = BedUtils.readBedRegions(bedFileName);
		mergeBedRegions(bedRegionMap, range);
		BedUtils.writeBedRegions(bedRegionMap.values()
				.stream()
				.flatMap(Collection::stream)
				.sorted(BedInterval::compareTo)
				.collect(Collectors.toList()), bedFileName + ".merged_range" + range);
	}


	private static void mergeBedRegions(Map<String, List<BedInterval>> bedRegionMap, int range) {
		for (Map.Entry<String, List<BedInterval>> bedRegionListEntry : bedRegionMap.entrySet()) {
			bedRegionListEntry.getValue().sort(BedInterval::compareTo);
			Iterator<BedInterval> iter = bedRegionListEntry.getValue().iterator();
			BedInterval curr = iter.next();
			while (iter.hasNext()) {
				BedInterval next = iter.next();
				if (next.getStart() - curr.getEnd() <= range) {
					curr.setEnd(next.getEnd() < curr.getEnd() ? curr.getEnd() : next.getEnd());
					curr.setName(curr.getName() + "," + next.getName());
					iter.remove();
				} else {
					curr = next;
				}
			}
		}
	}
}


