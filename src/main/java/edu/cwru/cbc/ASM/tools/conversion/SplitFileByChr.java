package edu.cwru.cbc.ASM.tools.conversion;

import org.apache.commons.cli.*;

import java.io.*;
import java.util.logging.Logger;

/**
 * Created by ke on 2/19/14.
 * Split original downloaded whole file into chromosomes.
 */
public class SplitFileByChr {
	private static final Logger logger = Logger.getLogger(SplitFileByChr.class.getName());

	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input file").build());
		options.addOption(Option.builder("o").hasArg().desc("output path").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputFileName = cmd.getOptionValue("a");
		String outputFilePath = cmd.getOptionValue("b");
		split(inputFileName, outputFilePath);
	}

	public static void split(String inputFileName, String outputFilePath) throws IOException {
		File outputPath = new File(outputFilePath);
		if (!outputPath.exists()) {
			if (!outputPath.mkdir()) {
				throw new RuntimeException("failed to create folder:" + outputFilePath);
			}
		}
		BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFileName));
		BufferedWriter bufferedWriter = null;
		String line, chr = "";
		String[] items;
		while ((line = bufferedReader.readLine()) != null) {
			if (line.startsWith("chrom") || line.startsWith("chr")) {
				continue;
			}
			items = line.split("\t");
			if (!chr.equals(items[0])) {
				if (bufferedWriter != null) {
					logger.info(chr + "\tfinished");
					bufferedWriter.close();
				}
				chr = items[0];
				bufferedWriter = new BufferedWriter(new FileWriter(String.format("%s/chr%s", outputFilePath, chr)));
				bufferedWriter.write(String.format("%s\n", line));
			} else {
				assert bufferedWriter != null;
				bufferedWriter.write(String.format("%s\n", line));
			}
		}
		if (bufferedWriter != null) {
			bufferedWriter.close();
		}
		bufferedReader.close();
	}
}
