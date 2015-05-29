package edu.cwru.cbc.ASM.intersect;

import com.google.common.base.Strings;
import edu.cwru.cbc.ASM.commons.GenomicInterval;
import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegionsPgm {
	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption("a", true, "bed file A");
		options.addOption("b", true, "bed file B");
		options.addOption("o", true, "output file name");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String aFileName = cmd.getOptionValue("a");
		String bFileName = cmd.getOptionValue("b");
		String outputFileName = cmd.getOptionValue("o");

		execution(aFileName, bFileName, outputFileName);
	}

	private static void execution(String aFileName, String bFileName, String outputFileName) throws IOException {
		Collection<GenomicInterval> regionsA = BedUtils.readSingleChromBedRegions(aFileName, true);
		Collection<GenomicInterval> regionsB = BedUtils.readSingleChromBedRegions(bFileName, true);
		Collection<GenomicInterval> intersections = BedUtils.intersect(regionsA, regionsB);
		BedUtils.writeBedRegions(intersections, outputFileName);
		BedUtils.writeBedWithIntersection(regionsA,
				new File(aFileName).getAbsolutePath() + "_intersected.bed");
		BedUtils.writeBedWithIntersection(regionsB,
				new File(bFileName).getAbsolutePath() + "_intersected.bed");

		long intersectedCount_regionA = regionsA.stream().filter(GenomicInterval::isIntersected).count();
		long nonIntersectedCount_regionA = regionsA.stream().filter(r -> !r.isIntersected()).count();
		long intersectedCount_regionB = regionsB.stream().filter(GenomicInterval::isIntersected).count();
		long nonIntersectedCount_regionB = regionsB.stream().filter(r -> !r.isIntersected()).count();
		int intersectionLength = intersections.stream().mapToInt(GenomicInterval::length).sum();
		int totalLengthRegionA = regionsA.stream().mapToInt(GenomicInterval::length).sum();
		int totalLengthRegionB = regionsB.stream().mapToInt(GenomicInterval::length).sum();

		System.out.println(Strings.padStart("", 80, '*'));
		System.out.printf("A is %s\n", aFileName);
		System.out.printf("B is %s\n", bFileName);

		System.out.println(Strings.padStart("", 80, '*'));
		System.out.printf("#intersection:%d\n", intersections.size());
		System.out.printf("total length of intersection:%d\n", intersectionLength);

		System.out.println(Strings.padStart("", 80, '*'));
		System.out.printf("total #regions in A:%d\n", regionsA.size());
		System.out.printf("total length of A:%d\n", totalLengthRegionA);
		System.out.printf("#covered regions in A:%d(%.4f)\n", intersectedCount_regionA,
				(double) intersectedCount_regionA / regionsA.size());
		System.out.printf("#non-covered regions in A:%d(%.4f)\n", nonIntersectedCount_regionA,
				(double) nonIntersectedCount_regionA / regionsA.size());
		System.out.printf("percent of intersected length of A:%.4f\n", (double) intersectionLength / totalLengthRegionA);
		assert intersectedCount_regionA + nonIntersectedCount_regionA == regionsA.size();

		System.out.println(Strings.padStart("", 80, '*'));
		System.out.printf("total #regions in B:%d\n", regionsB.size());
		System.out.printf("total length of B:%d\n", totalLengthRegionB);
		System.out.printf("#covered regions in B:%d(%.4f)\n", intersectedCount_regionB,
				(double) intersectedCount_regionB / regionsB.size());
		System.out.printf("#non-covered regions in B:%d(%.4f)\n", nonIntersectedCount_regionB,
				(double) nonIntersectedCount_regionB / regionsB.size());
		System.out.printf("percent of intersected length of B:%.4f\n", (double) intersectionLength / totalLengthRegionB);
		assert intersectedCount_regionB + nonIntersectedCount_regionB == regionsB.size();

		System.out.println(Strings.padStart("", 80, '*'));
	}
}
