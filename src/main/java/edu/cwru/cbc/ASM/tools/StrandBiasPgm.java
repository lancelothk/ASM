package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.io.GroupedReadsLineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import edu.cwru.cbc.ASM.detect.FisherExactTest;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by kehu on 7/6/16.
 * Program to calculate strand bias for reads from two alleles.
 */
public class StrandBiasPgm {
	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().required().desc("input grouped read file").build());
		options.addOption(Option.builder("n").desc("don't print header").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String groupedReadFile = cmd.getOptionValue("i");

		Pair<String, Pair<List<MappedRead>, List<MappedRead>>> result = Files.readLines(new File(groupedReadFile),
				Charsets.UTF_8, new GroupedReadsLineProcessor());
		List<MappedRead> group1 = result.getRight().getLeft();
		List<MappedRead> group2 = result.getRight().getRight();
		Pair<Double, Double> result1 = getStrandCounts(group1);
		Pair<Double, Double> result2 = getStrandCounts(group2);
		double a = result1.getLeft(), c = result1.getRight(), b = result2.getLeft(), d = result2.getRight();

		if (!cmd.hasOption("n")) {
			System.out.println("file\tStrandBias\tGATK SB\tFisher SB");
		}
		double sb1 = calculateStrandBias(a, b, c, d);
		double gatk1 = calculateGATKStrandBias(a, b, c, d);
		double fisher1 = calculateFisherStrandBias(a, b, c, d);
		System.out.printf("%s\t%f\t%f\t%f\t", groupedReadFile, sb1, gatk1, fisher1);

		a = result2.getLeft();
		c = result2.getRight();
		b = result1.getLeft();
		d = result1.getRight();
		double sb2 = calculateStrandBias(a, b, c, d);
		double gatk2 = calculateGATKStrandBias(a, b, c, d);
		double fisher2 = calculateFisherStrandBias(a, b, c, d);
		System.out.printf("%f\t%f\t%f\t", sb2, gatk2, fisher2);

		System.out.printf("%f\t%f\t%f\n", Math.max(sb1, sb2), Math.max(gatk1, gatk2), Math.max(fisher1, fisher2));
	}

	static double calculateStrandBias(double a, double b, double c, double d) {
		double ab = a + b, ac = a + c, ad = a + d, bc = b + c, bd = b + d, cd = c + d, abcd = a + b + c + d;
		return Math.abs(b / ab - d / cd) / (bd / (abcd));
	}

	static double calculateGATKStrandBias(double a, double b, double c, double d) {
		double ab = a + b, ac = a + c, ad = a + d, bc = b + c, bd = b + d, cd = c + d, abcd = a + b + c + d;
		return Math.max((b * c / ab / cd) / (ac / abcd), (d * a / cd / ab) / (ac / abcd));
	}

	static double calculateFisherStrandBias(double a, double b, double c, double d) {
		return 1 - FisherExactTest.fishersExactTest((int) a, (int) b, (int) c, (int) d)[0];
	}

	private static Pair<Double, Double> getStrandCounts(List<MappedRead> group1) {
		double g1_plus = 0, g1_minus = 0;
		for (MappedRead mappedRead : group1) {
			if (mappedRead.getStrand() == '+') {
				g1_plus++;
			} else if (mappedRead.getStrand() == '-') {
				g1_minus++;
			} else {
				throw new RuntimeException("unknown strand type!");
			}
		}
		return new ImmutablePair<>(g1_plus, g1_minus);
	}
}
