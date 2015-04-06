package edu.cwru.cbc.ASM.simulation;

import edu.cwru.cbc.ASM.commons.Bed.BedUtils;
import edu.cwru.cbc.ASM.commons.GenomicInterval;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegions {
	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption("a", true, "bed file A");
		options.addOption("b", true, "bed file B");
		options.addOption("o", true, "output path");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String aFileName = cmd.getOptionValue("a");
		String bFileName = cmd.getOptionValue("b");
		String outputPathName = cmd.getOptionValue("o");

		execution(aFileName, bFileName, outputPathName);
	}

	private static void execution(String aFileName, String bFileName, String outputPathName) throws IOException {
		Collection<GenomicInterval> regionsA = BedUtils.readSingleChromBedRegions(aFileName, true);
		Collection<GenomicInterval> regionsB = BedUtils.readSingleChromBedRegions(bFileName, true);
		Collection<GenomicInterval> intersections = BedUtils.intersect(regionsA, regionsB);
		BedUtils.writeBedRegions(intersections, outputPathName + "intersection.bed");
		BedUtils.writeBedWithIntersection(regionsA,
				outputPathName + new File(aFileName).getName() + "_intersected.bed");
		BedUtils.writeBedWithIntersection(regionsB,
				outputPathName + new File(bFileName).getName() + "_intersected.bed");

		double intersectedCount_regionA = regionsA.stream().filter(GenomicInterval::isIntersected).count();
		double nonIntersectedCount_regionA = regionsA.stream().filter(r -> !r.isIntersected()).count();
		double intersectedCount_regionB = regionsB.stream().filter(GenomicInterval::isIntersected).count();
		double nonIntersectedCount_regionB = regionsB.stream().filter(r -> !r.isIntersected()).count();
		int intersectionLength = intersections.stream().mapToInt(GenomicInterval::length).sum();
		int nonIntersectionLength_regionA =
				regionsA.stream().mapToInt(GenomicInterval::length).sum() - intersectionLength;
		int nonIntersectionLength_regionB =
				regionsB.stream().mapToInt(GenomicInterval::length).sum() - intersectionLength;

		System.out.printf("#intersection:%d\n", intersections.size());
		System.out.printf("non-intersected length of A:%d\n", nonIntersectionLength_regionA);
		System.out.printf("non-intersected length of B:%d\n", nonIntersectionLength_regionB);
		System.out.printf("total length of intersection:%d\n", intersectionLength);
		System.out.printf("#covered regions in A:%.4f(%.4f)\n", intersectedCount_regionA,
				intersectedCount_regionA / regionsA.size());
		System.out.printf("#not covered regions in A:%.4f(%.4f)\n", nonIntersectedCount_regionA,
				nonIntersectedCount_regionA / regionsA.size());
		System.out.printf("#regions in A:%d\n", regionsA.size());
		System.out.printf("#covered regions in B:%.4f(%.4f)\n", intersectedCount_regionB,
				intersectedCount_regionB / regionsB.size());
		System.out.printf("#not covered regions in B:%.4f(%.4f)\n", nonIntersectedCount_regionB,
				nonIntersectedCount_regionB / regionsB.size());
		System.out.printf("#regions in B:%d\n", regionsB.size());
	}

	public static void calcOverlappedLenth(String targetFileName, String actualResultFileName) throws IOException {
		int overlappedLength = 0, tpLength = 0;
		List<GenomicInterval> srcRegions = BedUtils.readSingleChromBedRegions(targetFileName);
		List<GenomicInterval> dstRegions = BedUtils.readSingleChromBedRegions(actualResultFileName);
		if (!srcRegions.get(0).getChr().equals(dstRegions.get(0).getChr())) {
			throw new RuntimeException("two Bed file should contain regions from same chromosome!");
		}
		for (GenomicInterval targetRegion : srcRegions) {
			for (GenomicInterval resultRegion : dstRegions) {
				if (targetRegion.getStart() <= resultRegion.getEnd() &&
						targetRegion.getEnd() >= resultRegion.getStart()) {
					targetRegion.addIntersectedRegion(resultRegion);
					resultRegion.addIntersectedRegion(targetRegion);
					int right = Math.min(targetRegion.getEnd(), resultRegion.getEnd());
					int left = Math.max(targetRegion.getStart(), resultRegion.getStart());
					int length = right - left + 1;
					System.out.println(resultRegion.getStart() + "\t" + resultRegion.getEnd() +
							"\toverlapped length:\t" + length + "\tnonoverlapped length:\t" +
							(resultRegion.length() - length));
					tpLength += (resultRegion.length() - length);
					overlappedLength += length;
				}
			}
		}
		System.out.println("overlapped length:\t" + overlappedLength);
		System.out.println("nonOverlapped length:\t" + (tpLength));
	}
}
