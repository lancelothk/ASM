package edu.cwru.cbc.ASM.tools.Simulation;

import edu.cwru.cbc.ASM.commons.DataType.GenomicRegion;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.readBedRegions;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegions {
	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption("t", true, "target file");
		options.addOption("r", true, "actual result file");
		options.addOption("e", true, "expected result file");
		options.addOption("o", true, "output file");
		options.addOption("a", true, "alpha");
		options.addOption("b", true, "beta");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String targetFileName = cmd.getOptionValue("t");
		String actualResultFileName = cmd.getOptionValue("r");
		String expectedResultFileName  = cmd.getOptionValue("e");
		String outputFileName = cmd.getOptionValue("o");
		double alpha  = Double.valueOf(cmd.getOptionValue("a"));
		double beta = Double.valueOf(cmd.getOptionValue("b"));

		execution(targetFileName, actualResultFileName, expectedResultFileName, outputFileName, alpha, beta);
	}

	public static void execution(String targetFileName, String actualResultFileName, String expectedResultFileName,
								 String outputFileName, double alpha, double beta) throws IOException {
		System.out.printf("%.1f-%.1f\n", alpha, beta);
		twoWayIntersection(targetFileName, actualResultFileName, expectedResultFileName);
		fourWayIntersection(expectedResultFileName, actualResultFileName, outputFileName);
		System.out.println();
	}

	private static void twoWayIntersection(String targetFileName, String actualResultFileName,
										   String expectedResultFileName) throws IOException {
		List<GenomicRegion> targetRegions = readBedRegions(targetFileName);
		List<GenomicRegion> resultRegions = readBedRegions(actualResultFileName);

		BufferedWriter writer = new BufferedWriter(new FileWriter(expectedResultFileName));

		for (GenomicRegion targetRegion : targetRegions) {
			for (GenomicRegion resultRegion : resultRegions) {
				if (targetRegion.getStart() <= resultRegion.getEnd() &&
						targetRegion.getEnd() >= resultRegion.getStart()) {
					targetRegion.setIntersected(true);
					resultRegion.setIntersected(true);
					writer.write(String.format("%s\t%s\t+\n", resultRegion.toBedString(), targetRegion.toBedString()));
				}
			}
		}

		for (GenomicRegion targetRegion : targetRegions) {
			if (!targetRegion.isIntersected()) {
				writer.write(targetRegion.toBedString() + "\t*" + "\n");
			}
		}

		for (GenomicRegion resultRegion : resultRegions) {
			if (!resultRegion.isIntersected()) {
				writer.write(resultRegion.toBedString() + "\t-" + "\n");
			}
		}

		writer.close();
	}

	private static void fourWayIntersection(String expectedResultFileName, String actualResultFileName,
											String outputFileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
		int tp = 0, fp = 0, fn = 0, tn = 0;
		List<GenomicRegion> expectedRegions = readBedRegions(expectedResultFileName, true);
		List<GenomicRegion> actualRegions = readBedRegions(actualResultFileName, true);
		for (GenomicRegion expectedRegion : expectedRegions) {
			for (GenomicRegion actualRegion : actualRegions) {
				if (expectedRegion.getStart() == actualRegion.getStart() &&
						expectedRegion.getEnd() == actualRegion.getEnd()) {
					if (expectedRegion.isPositive() && actualRegion.isPositive()) {
						tp++;
						writer.write(expectedRegion.toBedString() + "\ttp\n");
					}
					if (expectedRegion.isPositive() && !actualRegion.isPositive()) {
						fn++;
						writer.write(expectedRegion.toBedString() + "\tfn\n");
					}
					if (!expectedRegion.isPositive() && actualRegion.isPositive()) {
						fp++;
						writer.write(expectedRegion.toBedString() + "\tfp\n");
					}
					if (!expectedRegion.isPositive() && !actualRegion.isPositive()) {
						tn++;
						writer.write(expectedRegion.toBedString() + "\ttn\n");
					}
				}
			}
		}
		writer.close();
		System.out.printf("truePositive:%d\ttrueNegative:%d\n", tp, tn);
		System.out.printf("falseNegative:%d\tfalsePositive:%d\n", fn, fp);
		System.out.printf("accuracy:%.4f\n", (tp + tn) / (double) (tp + tn + fn + fp));
		System.out.printf("Precision:%.4f\n", tp / (double) (tp + fp));
		System.out.printf("Sensitivity:%.4f\n", tp / (double) (tp + fn));
		System.out.printf("F1 score:%.4f\n", 2 * tp / (double) (2 * tp + fp + fn));
		System.out.printf("MCC:%.4f\n", (tp * tn - fp * fn) / (double) ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)));
	}
}
