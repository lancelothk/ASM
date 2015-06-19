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
		options.addOption("i", true, "Input intervals folder or interval file name");
		options.addOption("mic", true, "Minimum interval cpg number");
		options.addOption("mcc", true, "Minimum adjacent CpG coverage");
		options.addOption("f", true, "FDR threshold");
		options.addOption("t", false, "Thread number to execute the program.");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String inputPath = cmd.getOptionValue("i");
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		int min_cpg_coverage = Integer.valueOf(cmd.getOptionValue("mcc"));
		double FDR_threshold = Double.valueOf(cmd.getOptionValue("f"));
		int threadNumber = Integer.valueOf(cmd.getOptionValue("t", "6"));
		execute(inputPath, threadNumber, min_interval_cpg, min_cpg_coverage, FDR_threshold);
		System.out.println(System.currentTimeMillis() - start + "ms");
	}

	private static void execute(String inputPath, int threadNumber, int min_interval_cpg, int min_cpg_coverage,
	                            double FDR_threshold) throws ExecutionException, InterruptedException, IOException {
		// initialize IntervalDetectionSummary format
		IntervalDetectionSummary.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>().add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#refCpG", "%d"))
						.add(new ImmutablePair<>("#clusterCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("errorProb", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("group1", "%d"))
						.add(new ImmutablePair<>("group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.add(new ImmutablePair<>("label", "%s"))
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
									new Detection(file, min_interval_cpg, min_cpg_coverage));
							futureList.add(future);
						}
					} catch (Exception e) {
						throw new RuntimeException("Problem File name: " + file.getAbsolutePath() + "\n", e);
					}
				}
			}
		} else {
			try {
				Future<IntervalDetectionSummary> future = executor.submit(
						new Detection(inputFile, min_interval_cpg, min_cpg_coverage));
				futureList.add(future);
			} catch (Exception e) {
				throw new RuntimeException("Problem File name:" + inputFile.getAbsolutePath(), e);
			}
		}
		for (Future<IntervalDetectionSummary> intervalDetectionSummaryFuture : futureList) {
			resultList.add(intervalDetectionSummaryFuture.get());
		}
		executor.shutdown();

		double regionP_threshold = FDRControl.getBHYFDRCutoff(
				resultList.stream()
						.map(IntervalDetectionSummary::getRegionP)
						.filter(p -> p <= 1 && p >= 0)
						.collect(Collectors.toList()),
				FDR_threshold);
		System.out.println("regionP threshold calculated by FDR control:\t" + regionP_threshold);
		writeDetectionSummary(inputPath, resultList, regionP_threshold);
	}

	private static void writeDetectionSummary(String outputPath, List<IntervalDetectionSummary> resultList,
	                                          double region_threshold) throws IOException {
		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(outputPath + "/detection.summary"));
		summaryWriter.write(IntervalDetectionSummary.getHeadLine());
		for (IntervalDetectionSummary result : resultList) {
			if (result.getRegionP() <= region_threshold) {
				summaryWriter.write(result.getSummaryString(region_threshold));
			}
		}
		summaryWriter.close();
	}
}
