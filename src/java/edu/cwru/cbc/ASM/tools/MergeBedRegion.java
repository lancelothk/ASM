package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;

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

	public static void main(String[] args) throws IOException {
		int range = 1000;
		String bedFileName = "/home/lancelothk/result/h1.detection.summary.combined.chr.bed";
		Map<String, List<BedInterval>> bedRegionMap = BedUtils.readBedRegions(bedFileName);
		mergeBedRegions(bedRegionMap, range);
		BedUtils.writeBedRegions(bedRegionMap.values()
				.stream()
				.flatMap(Collection::stream)
				.sorted(BedInterval::compareTo)
				.collect(Collectors.toList()), bedFileName + ".merged");
	}

	private static void mergeBedRegions(Map<String, List<BedInterval>> bedRegionMap, int range) {
		for (Map.Entry<String, List<BedInterval>> bedRegionListEntry : bedRegionMap.entrySet()) {
			bedRegionListEntry.getValue().sort(BedInterval::compareTo);
			Iterator<BedInterval> iter = bedRegionListEntry.getValue().iterator();
			BedInterval curr = iter.next();
			while (iter.hasNext()) {
				BedInterval next = iter.next();
				if (next.getStart() - curr.getEnd() <= range) {
					curr.setEnd(next.getEnd());
//					curr.setName(curr.getName() + "," + next.getName());
					iter.remove();
				} else {
					curr = next;
				}
			}
		}
	}
}


