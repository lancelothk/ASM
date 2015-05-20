package edu.cwru.cbc.ASM.intersect;

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
}
