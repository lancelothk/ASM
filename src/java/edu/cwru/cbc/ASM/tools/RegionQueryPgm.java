package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by kehu on 9/16/15.
 */
public class RegionQueryPgm {
	private static final Splitter tabSplitter = Splitter.on("\t");

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input File").required().build());
		options.addOption(Option.builder("c").hasArg().desc("chromosome of input").required().build());
		options.addOption(Option.builder("s").hasArg().desc("start").required().build());
		options.addOption(Option.builder("e").hasArg().desc("end").required().build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String inputFileName = cmd.getOptionValue("i");
		String inputChr = cmd.getOptionValue("c");
		int startPos = Integer.parseInt(cmd.getOptionValue("s"));
		int endPos = Integer.parseInt(cmd.getOptionValue("e"));


		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor() {
			private int CpGCount = 0;
			private double sumCoverage = 0;
			private double sumMethylLevel = 0;
			private double weightedSumMethylLevel = 0;

			@Override
			public boolean processLine(String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				String chr = itemList.get(0);
				int pos = Integer.parseInt(itemList.get(1));
				int coverage = Integer.parseInt(itemList.get(2));
				double methylLevel = Double.parseDouble(itemList.get(4));
				if (inputChr.equals(chr) && pos >= startPos && pos <= endPos && coverage > 0) {
					CpGCount++;
					sumCoverage += coverage;
					sumMethylLevel += methylLevel;
					weightedSumMethylLevel += coverage * methylLevel;
				}
				return true;
			}

			@Override
			public Object getResult() {
				System.out.printf("%s\t%d\t%d\t%d\t%f\t%f\t%f\n", inputChr, startPos, endPos, CpGCount,
						sumCoverage / CpGCount, sumMethylLevel / CpGCount, weightedSumMethylLevel / sumCoverage);
				return null;
			}
		});
	}
}