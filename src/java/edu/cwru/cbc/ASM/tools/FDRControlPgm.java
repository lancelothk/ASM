package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.detect.FDRControl;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.File;
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
		options.addOption(Option.builder("i").hasArg().desc("input grouped read file").build());
		options.addOption(Option.builder("c").hasArg().desc("column of p values").build());
		options.addOption(Option.builder("f").hasArg().desc("fdr threshold").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String inputFile = cmd.getOptionValue("i");
		double fdr = Double.parseDouble(cmd.getOptionValue("f"));
		int column = Integer.parseInt(cmd.getOptionValue("c")) - 1;

		double correctedThreshold = Files.readLines(new File(inputFile), Charsets.UTF_8, new LineProcessor<Double>() {
			private List<Double> pvalueList = new ArrayList<>();

			@Override
			public boolean processLine(@Nonnull String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				if (line.startsWith("ref:") || line.equals("") || line.startsWith("chr\t")) {
					return true;
				} else {
					pvalueList.add(Double.parseDouble(itemList.get(column)));
					return true;
				}
			}

			@Override
			public Double getResult() {
				return FDRControl.getBHYFDRCutoff(pvalueList, fdr);
			}
		});

		System.out.println("corrected threshold:" + correctedThreshold);
	}
}
