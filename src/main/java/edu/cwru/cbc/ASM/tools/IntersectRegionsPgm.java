package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Strings;
import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegionsPgm {
	public static void main(String[] args) throws IOException, ParseException, ExecutionException,
			InterruptedException {
		Options options = new Options();
		options.addOption(Option.builder("a").hasArg().desc("bed file A").required().build());
		options.addOption(Option.builder("b").hasArg().desc("bed file B").required().build());
		options.addOption(Option.builder("o")
				.hasArg()
				.desc("output file name prefix(without .bed)")
				.required()
				.build());
		options.addOption(Option.builder("e")
				.hasArg()
				.desc("extra output except intersection file. 'a':annotated A; 'b':annotated B;'c':both annotated A&B ")
				.required(false)
				.build());
		options.addOption(Option.builder("t").hasArg().required(false).desc("number of threads to use").build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String aFileName = cmd.getOptionValue("a");
		String bFileName = cmd.getOptionValue("b");
		String outputFileName = cmd.getOptionValue("o");
		int threadNumber = Integer.parseInt(cmd.getOptionValue("t", "1"));
		String extraOption = "";
		if (cmd.hasOption("e")) {
			extraOption = cmd.getOptionValue("e");
			if (!extraOption.equals("a") && !extraOption.equals("b") && !extraOption.equals("c")) {
				System.err.println(
						"unknown parameter \"" + extraOption + "\" for -e.\n'a':annotated A; 'b':annotated B;'c':both annotated A&B");
			}
		}
		execution(aFileName, bFileName, outputFileName, extraOption, threadNumber);
	}

	private static void execution(String aFileName, String bFileName, String outputFileName, String extraOption,
	                              int threadNumber) throws
			IOException, ExecutionException, InterruptedException {
		Map<String, List<BedInterval>> regionsMapA = BedUtils.readBedRegions(aFileName);
		Map<String, List<BedInterval>> regionsMapB = BedUtils.readBedRegions(bFileName);
		List<BedInterval> intersections = new ArrayList<>();
		if (threadNumber == 1) {
			regionsMapA.keySet().stream().filter(regionsMapB::containsKey).forEach(keyA ->
					intersections.addAll(BedUtils.intersect(regionsMapA.get(keyA), regionsMapB.get(keyA))));
		} else {
			ForkJoinPool forkJoinPool = new ForkJoinPool(threadNumber);
			forkJoinPool.submit(
					() -> regionsMapA.keySet().stream().parallel().filter(regionsMapB::containsKey).forEach(keyA ->
							intersections.addAll(BedUtils.intersect(regionsMapA.get(keyA), regionsMapB.get(keyA))))
			).get();
		}
		BedUtils.writeBedRegions(intersections, outputFileName + ".bed");
		List<BedInterval> regionsListA = regionsMapA.values().stream().flatMap(Collection::stream).sorted(
				BedInterval::compareTo).collect(Collectors.toList());
		List<BedInterval> regionsListB = regionsMapB.values().stream().flatMap(Collection::stream).sorted(
				BedInterval::compareTo).collect(Collectors.toList());
		switch (extraOption) {
			case "a":
				BedUtils.writeBedWithIntersection(regionsListA, outputFileName + "_A.bed");
				break;
			case "b":
				BedUtils.writeBedWithIntersection(regionsListB, outputFileName + "_B.bed");
				break;
			case "c":
				BedUtils.writeBedWithIntersection(regionsListA, outputFileName + "_A.bed");
				BedUtils.writeBedWithIntersection(regionsListB, outputFileName + "_B.bed");
				break;
			default:
		}
		printSummary(aFileName, bFileName, regionsListA, regionsListB, intersections);
	}

	private static void printSummary(String aFileName, String bFileName, Collection<BedInterval> regionsA, Collection<BedInterval> regionsB, List<BedInterval> intersections) {
		long intersectedCount_regionA = regionsA.stream().filter(BedInterval::isIntersected).count();
		long nonIntersectedCount_regionA = regionsA.stream().filter(r -> !r.isIntersected()).count();
		long intersectedCount_regionB = regionsB.stream().filter(BedInterval::isIntersected).count();
		long nonIntersectedCount_regionB = regionsB.stream().filter(r -> !r.isIntersected()).count();
		long intersectionLength = intersections.stream().mapToLong(BedInterval::length).sum();
		long totalLengthRegionA = regionsA.stream().mapToLong(BedInterval::length).sum();
		long totalLengthRegionB = regionsB.stream().mapToLong(BedInterval::length).sum();

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
		System.out.printf("percent of intersected length of A:%.4f\n",
				(double) intersectionLength / totalLengthRegionA);
		assert intersectedCount_regionA + nonIntersectedCount_regionA == regionsA.size();

		System.out.println(Strings.padStart("", 80, '*'));
		System.out.printf("total #regions in B:%d\n", regionsB.size());
		System.out.printf("total length of B:%d\n", totalLengthRegionB);
		System.out.printf("#covered regions in B:%d(%.4f)\n", intersectedCount_regionB,
				(double) intersectedCount_regionB / regionsB.size());
		System.out.printf("#non-covered regions in B:%d(%.4f)\n", nonIntersectedCount_regionB,
				(double) nonIntersectedCount_regionB / regionsB.size());
		System.out.printf("percent of intersected length of B:%.4f\n",
				(double) intersectionLength / totalLengthRegionB);
		assert intersectedCount_regionB + nonIntersectedCount_regionB == regionsB.size();

		System.out.println(Strings.padStart("", 80, '*'));
	}
}
