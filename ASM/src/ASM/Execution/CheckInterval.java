package ASM.Execution;

import ASM.DataType.ChrCoverageSummary;
import ASM.Utils.IntervalChekingLineProcessor;
import com.google.common.io.CharStreams;
import org.apache.commons.cli.*;

import java.io.*;
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
			System.out.println("\nArgs corrrect?y/n");
			Scanner check = new Scanner(System.in);
			String answer = check.nextLine();
			if (answer.equals("y")) {
				checkInterval(chrSize, chr, referenceFileName, mappedReadFileName, outputFileName);
			}
		} catch (ParseException e) {
			e.printStackTrace();
		}
	}

	public static void checkInterval(int chrBitSize, String chr, String referenceFileName, String mappedReadFileName, String outputFileName) throws IOException {
		ChrCoverageSummary chrCoverageSummary = CharStreams.readLines(new BufferedReader((new FileReader(mappedReadFileName))), new IntervalChekingLineProcessor(chrBitSize, referenceFileName));
		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		chrCoverageSummary.writeIntervalCoverageSummary(chr, bufferedWriter);
		bufferedWriter.close();
	}
}
