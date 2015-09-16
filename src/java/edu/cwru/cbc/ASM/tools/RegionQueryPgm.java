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
		options.addOption("i", true, "input File");
		options.addOption("s", true, "start pos");
		options.addOption("e", true, "end pos(inclusive)");
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String inputFileName = cmd.getOptionValue("i");
		int startPos = Integer.parseInt(cmd.getOptionValue("s"));
		int endPos = Integer.parseInt(cmd.getOptionValue("e"));


		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor() {
			private int CpGCount = 0;
			private double avgCoverage = 0;
			private double avgmethylLevel = 0;
			private double weightedAvgMethylLevel = 0;

			@Override
			public boolean processLine(String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				int pos = Integer.parseInt(itemList.get(0));
				int coverage = Integer.parseInt(itemList.get(1));
				double methylLevel = Double.parseDouble(itemList.get(3));
				if (pos >= startPos && pos <= endPos && coverage > 0) {
					CpGCount++;
					avgCoverage += coverage;
					avgmethylLevel += methylLevel;
					weightedAvgMethylLevel += coverage * methylLevel;
				}
				return true;
			}

			@Override
			public Object getResult() {
				System.out.println("cpg count:" + CpGCount);
				System.out.println("avg coverage:" + avgCoverage / CpGCount);
				System.out.println("avg methyl level:" + avgmethylLevel / CpGCount);
				System.out.println("weighted avg methyl level:" + weightedAvgMethylLevel / avgCoverage);
				return null;
			}
		});
	}
}