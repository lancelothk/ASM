package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.detect.dataType.IntervalDetectionSummaryFormatter;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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
		options.addOption(Option.builder("p").hasArg().desc("Time of random permutation").required().build());
		options.addOption(Option.builder("t").hasArg().desc("Thread number to call the program").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputPath = cmd.getOptionValue("i");
		int min_interval_cpg = Integer.valueOf(cmd.getOptionValue("mic"));
		int permTime = Integer.valueOf(cmd.getOptionValue("p"));
		int threadNumber = Integer.valueOf(cmd.getOptionValue("t", "1"));
		execute(inputPath, threadNumber, min_interval_cpg, permTime);
		System.out.println(System.currentTimeMillis() - start + "ms");
	}

	private static void execute(String inputPath, int threadNumber, int min_interval_cpg, int permTime) throws
			ExecutionException, InterruptedException, IOException {
		// initialize IntervalDetectionSummaryFormatter format
		IntervalDetectionSummaryFormatter.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>()
						.add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("originalRegion", "%s"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#clusterRefCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("clusterIndex", "%f"))
						.add(new ImmutablePair<>("#group1", "%d"))
						.add(new ImmutablePair<>("#group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.build());

		File inputFile = new File(inputPath);
		List<String> resultList = new ArrayList<>();
		ExecutorService executor = Executors.newFixedThreadPool(threadNumber);
		List<Future<String>> futureList = new ArrayList<>();
		if (inputFile.isDirectory()) {
			File[] files = inputFile.listFiles(
					pathname -> pathname.isFile() && pathname.getName().endsWith(Constant.MAPPEDREADS_EXTENSION));
			if (files == null) {
				throw new RuntimeException("Empty folder!");
			} else {
				Arrays.sort(files, (o1, o2) -> Long.compare(o2.length(), o1.length()));
				for (File file : files) {
					try {
						Future<String> future = executor.submit(
								new Detection(file, min_interval_cpg, permTime));
						futureList.add(future);
					} catch (Exception e) {
						throw new RuntimeException("Problem File name: " + file.getAbsolutePath() + "\n", e);
					}
				}
				for (Future<String> intervalDetectionSummaryFuture : futureList) {
					resultList.add(intervalDetectionSummaryFuture.get());
				}
				executor.shutdown();
			}
		} else {
			try {
				Future<String> future = executor.submit(
						new Detection(inputFile, min_interval_cpg, permTime));
				futureList.add(future);
			} catch (Exception e) {
				throw new RuntimeException("Problem File name:" + inputFile.getAbsolutePath(), e);
			}
			for (Future<String> intervalDetectionSummaryFuture : futureList) {
				resultList.add(intervalDetectionSummaryFuture.get());
			}
			executor.shutdown();
		}
		writeDetectionSummary(inputPath, resultList);
	}

	private static void writeDetectionSummary(String outputPath, List<String> resultList) throws IOException {
		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(outputPath + "/detection.summary"));
		summaryWriter.write(IntervalDetectionSummaryFormatter.getHeadLine() + "\n");
		for (String result : resultList) {
			if (!result.equals("")) {
				summaryWriter.write(result + "\n");
			}
		}
		summaryWriter.close();
	}
}
