package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.detect.FDRControl;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 10/26/15.
 * Program to calculate corrected FDR threshold for given p-values.
 */
public class FDRControlPgm {
	private static final Splitter tabSplitter = Splitter.on("\t");

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input summary file").build());
		options.addOption(Option.builder("o").hasArg().desc("output file with region pass FDR control").build());
		options.addOption(Option.builder("c").hasArg().desc("column of p values").build());
		options.addOption(Option.builder("f").hasArg().desc("fdr threshold").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String inputFile = cmd.getOptionValue("i");
		String outputFile = cmd.getOptionValue("o");
		double fdr = Double.parseDouble(cmd.getOptionValue("f"));
		int column = Integer.parseInt(cmd.getOptionValue("c")) - 1;

		Files.readLines(new File(inputFile), Charsets.UTF_8, new LineProcessor<Double>() {
			private List<Double> pvalueList = new ArrayList<>();
			private List<String> lines = new ArrayList<>();

			@Override
			public boolean processLine(@Nonnull String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				if (line.startsWith("ref:") || line.equals("") || line.startsWith("chr\t")) {
					return true;
				} else {
					pvalueList.add(Double.parseDouble(itemList.get(column)));
					lines.add(line);
					return true;
				}
			}

			@Override
			public Double getResult() {
				double correctedThreshold = FDRControl.getBHYFDRCutoff(pvalueList, fdr);

				try {
					BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
					for (int i = 0; i < lines.size(); i++) {
						if (pvalueList.get(i) <= correctedThreshold) {
							outputWriter.write(lines.get(i) + "\t");
						}
					}
					outputWriter.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}

				System.out.println("corrected threshold:" + correctedThreshold);
				return null;
			}
		});


	}
}
