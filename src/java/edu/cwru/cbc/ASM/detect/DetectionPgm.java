package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.detect.dataType.IntervalDetectionSummary;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

/**
 * Created by kehu on 11/12/14.
 * ASM Detection with whole read info.
 */
public class DetectionPgm {

	public static void main(
			String[] args) throws ParseException, IOException, ExecutionException, InterruptedException {
		long start = System.currentTimeMillis();

		Options options = new Options();
		options.addOption(
				Option.builder("i").hasArg().desc("Input intervals folder or interval file name").required().build());
		options.addOption(Option.builder("mic").hasArg().desc("Minimum interval CpG number").required().build());
		options.addOption(Option.builder("f").hasArg().desc("FDR threshold").required().build());
		options.addOption(
				Option.builder("p").hasArg().desc("p value threshold for fisher exact test").required().build());
		options.addOption(Option.builder("r").hasArg().desc("Time of random permutation").required().build());
		options.addOption(Option.builder("t").hasArg().desc("Thread number to call the program").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputPath = cmd.getOptionValue("i");
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		double FDR_threshold = Double.valueOf(cmd.getOptionValue("f"));
		double fisher_p = Double.valueOf(cmd.getOptionValue("p"));
		int permTime = Integer.valueOf(cmd.getOptionValue("r"));
		int threadNumber = Integer.valueOf(cmd.getOptionValue("t", "1"));
		execute(inputPath, threadNumber, min_interval_cpg, FDR_threshold, fisher_p, permTime);
		System.out.println(System.currentTimeMillis() - start + "ms");
	}

	private static void execute(String inputPath, int threadNumber, int min_interval_cpg, double FDR_threshold,
	                            double fisher_p, int permTime) throws ExecutionException, InterruptedException,
			IOException {
		// initialize IntervalDetectionSummary format
		IntervalDetectionSummary.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>()
						.add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("originalRegion", "%s"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#refCpG", "%d"))
						.add(new ImmutablePair<>("#clusterCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("randPIndex", "%d"))
						.add(new ImmutablePair<>("clusterIndex", "%f"))
						.add(new ImmutablePair<>("#group1", "%d"))
						.add(new ImmutablePair<>("#group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.add(new ImmutablePair<>("FDR_label", "%s"))
						.build());

		File inputFile = new File(inputPath);
		List<IntervalDetectionSummary> resultList = new ArrayList<>();
		ExecutorService executor = Executors.newFixedThreadPool(threadNumber);
		List<Future<IntervalDetectionSummary>> futureList = new ArrayList<>();
		if (inputFile.isDirectory()) {
			File[] files = inputFile.listFiles();
			if (files == null) {
				throw new RuntimeException("Empty folder!");
			} else {
				for (File file : files) {
					try {
						if (file.isFile() && file.getName().endsWith(Constant.MAPPEDREADS_EXTENSION)) {
							Future<IntervalDetectionSummary> future = executor.submit(
									new Detection(file, min_interval_cpg, fisher_p, permTime));
							futureList.add(future);
						}
					} catch (Exception e) {
						throw new RuntimeException("Problem File name: " + file.getAbsolutePath() + "\n", e);
					}
				}
				for (Future<IntervalDetectionSummary> intervalDetectionSummaryFuture : futureList) {
					resultList.add(intervalDetectionSummaryFuture.get());
				}
				executor.shutdown();

				resultList = resultList.stream()
						.filter(i -> i.getRegionP() <= 1 && i.getRegionP() >= 0 && !i.isRandom())
						.collect(Collectors.toList());
				Optional<Double> max = resultList.stream().map(IntervalDetectionSummary::getRegionP).max(
						Double::compareTo);
				double regionP_threshold = max.isPresent() ? max.get() : -1;
				if (FDR_threshold < 0) {
					System.out.println("skip FDR control");
				} else {
					regionP_threshold = FDRControl.getBHYFDRCutoff(
							resultList.stream().map(IntervalDetectionSummary::getRegionP)
									.collect(Collectors.toList()), FDR_threshold);
					System.out.println("regionP threshold calculated by FDR control:\t" + regionP_threshold);
				}
				writeDetectionSummary(inputPath, resultList, regionP_threshold);
			}
		} else {
			try {
				Future<IntervalDetectionSummary> future = executor.submit(
						new Detection(inputFile, min_interval_cpg, fisher_p, permTime));
				futureList.add(future);
			} catch (Exception e) {
				throw new RuntimeException("Problem File name:" + inputFile.getAbsolutePath(), e);
			}
			for (Future<IntervalDetectionSummary> intervalDetectionSummaryFuture : futureList) {
				resultList.add(intervalDetectionSummaryFuture.get());
			}
			System.out.println(resultList.get(0).getSummaryString(-1));
			executor.shutdown();
		}

	}

	private static void writeDetectionSummary(String outputPath, List<IntervalDetectionSummary> resultList,
	                                          double region_threshold) throws IOException {
		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(outputPath + "/detection.summary"));
		summaryWriter.write(IntervalDetectionSummary.getHeadLine());
		for (IntervalDetectionSummary result : resultList) {
			if (result.getRegionP() <= region_threshold) {
				summaryWriter.write("chr" + result.getSummaryString(region_threshold));
			}
		}
		summaryWriter.close();
	}
}
