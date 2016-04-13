package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Align reads to readable format with padding.
 */
public class ReadsVisualizationPgm {
	private static String ref;

	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption("i", true, "input file");

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputFileName = cmd.getOptionValue("i");
		ReadsVisualizationPgm.align(inputFileName);
	}

	private static void align(String fileName) throws IOException {
		File targetFile = new File(fileName);
		if (targetFile.isDirectory()) {
			File[] files = targetFile.listFiles();
			assert files != null;
			for (File file : files) {
				if (file.isFile() && !file.isHidden() && !file.getName().endsWith(".aligned")) {
					String outputFileName = file.getAbsolutePath() + ".aligned";
					List<MappedRead> readsList = readMappedReads(file);
					writeAlignedReads(readsList, ref, outputFileName);
				}
			}
		} else {
			String outputFileName = fileName + ".aligned";
			List<MappedRead> readsList = readMappedReads(targetFile);
			writeAlignedReads(readsList, ref, outputFileName);
		}
	}

	private static List<MappedRead> readMappedReads(File inputFile) throws IOException {
		ref = IOUtils.readRefFromIntervalReadsFile(inputFile);
		return Files.readLines(inputFile, Charsets.UTF_8, new MappedReadLineProcessor());
	}

	private static void writeAlignedReads(List<MappedRead> readsList, String ref, String outputFileName) throws
			IOException {
		// sort reads first
		readsList.sort(MappedRead::compareTo);
		// set initial position
		int initialPos = readsList.get(0).getStart();
		BufferedWriter bufferedWriter;
		bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		if (ref != null) {
			bufferedWriter.write(String.format("ref:\t%s\n", ref));
		}
		for (MappedRead read : readsList) {
			bufferedWriter.write(read.toVisualizationString(initialPos) + "\n");
		}
		bufferedWriter.close();
	}

	public static void writeAlignedReadsIntoGroups(List<List<MappedRead>> readGroups, String ref, String outputFileName) throws
			IOException {
		List<MappedRead> allReads = new ArrayList<>();
		readGroups.forEach(allReads::addAll);
		// sort reads first
		allReads.sort(MappedRead::compareTo);
		// set initial position
		int initialPos = allReads.get(0).getStart();
		BufferedWriter bufferedWriter;
		bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		if (ref != null) {
			bufferedWriter.write(String.format("ref:\t%s\n", ref));
		}
		for (List<MappedRead> readGroup : readGroups) {
			readGroup.sort(MappedRead::compareTo);
			for (MappedRead read : readGroup) {
				bufferedWriter.write(read.toVisualizationString(initialPos) + "\n");
			}
			bufferedWriter.write("\n");
		}
		bufferedWriter.close();
	}
}