package ASM.Execution;

import ASM.DataType.ChrCoverageSummary;
import ASM.DataType.GenomicInterval;
import ASM.DataType.IntervalMatchingLineProcessor;
import ASM.Utils.IntervalChekingLineProcessor;
import com.google.common.io.CharStreams;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.List;
import java.util.Scanner;

/**
 * Created by lancelothk on 2/17/14.
 * Check mapped reads interval in Chromosome
 */
public class CheckInterval {
//	public static final int CHR6SIZE = 170899992;

	public static void main(String[] args) throws IOException {
		// TODO make options required. Check options more carefully.
		Options options = new Options();
		options.addOption(new Option("c", true, "chromosome name, like chr6"));
		options.addOption(new Option("s", true, "chromosome size, a number"));
		options.addOption(new Option("r", true, "reference file"));
		options.addOption(new Option("i", true, "input file"));
		options.addOption(new Option("o", true, "output file"));
        options.addOption(new Option("y", false, "yes to execute with current parameter"));

		CommandLineParser parser = new GnuParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			String chr = cmd.getOptionValue("c");
			String chrSizeStr = cmd.getOptionValue("s");
			if (chrSizeStr == null) {
				throw new RuntimeException("chromosome size is null!");
			}
			int chrSize = Integer.parseInt(chrSizeStr);
			String referenceFileName = cmd.getOptionValue("r");
			String mappedReadFileName = cmd.getOptionValue("i");
			String outputFileName = cmd.getOptionValue("o");
			for (String arg : args) {
				System.out.print(arg + "\t");
			}
            if (cmd.hasOption("y")) {
                checkInterval(chrSize, chr, referenceFileName, mappedReadFileName, outputFileName);
            }else{
                System.out.println("\nArgs corrrect?y/n");
                Scanner check = new Scanner(System.in);
                String answer = check.nextLine();
                if (answer.equals("y")) {
                    checkInterval(chrSize, chr, referenceFileName, mappedReadFileName, outputFileName);
                }
            }
		} catch (ParseException e) {
			e.printStackTrace();
		}
	}

	public static void checkInterval(int chrBitSize, String chr, String referenceFileName, String mappedReadFileName, String outputFileName) throws IOException {
		ChrCoverageSummary chrCoverageSummary = CharStreams.readLines(
				new BufferedReader((new FileReader(mappedReadFileName))), new IntervalChekingLineProcessor(chrBitSize));
        List<GenomicInterval> intervalList = CharStreams.readLines(
				new BufferedReader((new FileReader(mappedReadFileName))),
				new IntervalMatchingLineProcessor(chrCoverageSummary));
		writeIntervalSummary(outputFileName, chr, referenceFileName, intervalList);
	}

	private static void writeIntervalSummary(String outputFileName, String chr, String referenceFileName, List<GenomicInterval> intervalList) throws IOException {
		String reference = readReference(referenceFileName);
		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		bufferedWriter.write("chr\tstart\tend\tlength\treadCount\tmaxCount\tavgCount\tCpGCount\n");
		for (GenomicInterval genomicInterval : intervalList) {
			bufferedWriter.write(String.format("%s\t%d\t%d\t%d\t%d\t%d\t%f\t%d\n", chr, genomicInterval.getStart(),
					genomicInterval.getEnd(), genomicInterval.getLength(), genomicInterval.getReadCount(),
					genomicInterval.getMaxCount(), genomicInterval.getAvgCount(),
					countCpG(reference, genomicInterval.getStart(), genomicInterval.getEnd())));
		}
		bufferedWriter.close();
	}

	private static String readReference(String fileName) throws IOException {
		List<String> lines = CharStreams.readLines(new BufferedReader(new FileReader(fileName)));
		StringBuilder referenceBuilder = new StringBuilder();
		lines.remove(0); // remove first line, which is chromosome name
		for (String line : lines) {
			referenceBuilder.append(line);
		}
		return referenceBuilder.toString();
	}

	private static int countCpG(String reference, int start, int end) {
		// start and end is 1-based.
		int count = 0;
		for (int i = start - 1; i <= end - 1; i++) {
			if (i + 1 < reference.length() && (reference.charAt(i) == 'c' || reference.charAt(
					i) == 'C') && (reference.charAt(i + 1) == 'g' || reference.charAt(
					i + 1) == 'G')) { // 'C' is last character in reference
				count++;
			}
		}
		return count;
	}
}
