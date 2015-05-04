package edu.cwru.cbc.ASM.conversion;

import edu.cwru.cbc.ASM.conversion.Conversion.MRToMappedRead;
import org.apache.commons.cli.*;

import java.io.IOException;

import static edu.cwru.cbc.ASM.conversion.Conversion.MappedReadToMR.conversion;

/**
 * Created by kehu on 5/4/15.
 * FormatConversion program entrance.
 */
public class FormatConversion {
	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption("m", true, "conversion mode. 1:MappedRead->MR;2.MR->MappedRead");
		options.addOption("i", true, "input MappedRead file");
		options.addOption("o", true, "output MR file");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String mode = cmd.getOptionValue("m");
		String inputFileName = cmd.getOptionValue("i");
		String outputFileName = cmd.getOptionValue("o");

		switch (mode) {
			case "1":
				conversion(inputFileName, outputFileName);
				break;
			case "2":
				MRToMappedRead.conversion(inputFileName, outputFileName);
				break;
			default:
				System.err.println("Unknown mode!");
		}

	}
}
