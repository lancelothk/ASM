package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;

/**
 * Created by lancelothk on 6/16/15.
 * Used for parsing experiment intersection result.
 */
public class ParseResult {
	public static void main(String[] args) throws IOException, ParseException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input result file").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String resultFileName = cmd.getOptionValue("i");
		String result = Files.readLines(new File(resultFileName), Charsets.UTF_8,
				new LineProcessor<String>() {
					private StringBuilder sb = new StringBuilder();

					@Override
					public boolean processLine(String line) throws IOException {
						if (line.startsWith("CGI") || line.startsWith("nonCGI")) {
							sb.append(line).append("-");
						} else if (line.startsWith("1")) {
							sb.append(line).append("\n");
						} else if (line.startsWith("0")) {
							sb.append(line).append("\t");
						} else if (line.startsWith("#covered regions in")) {
							sb.append(line.split(":")[1].split("\\(")[0]).append("\t");
							sb.append(line.split(":")[1].split("\\(")[1].replace(")", "")).append("\t");
						} else if (line.startsWith("total length of intersection")) {
							sb.append(line.split(":")[1]).append("\t");
						} else if (line.startsWith("percent of intersected length of")) {
							sb.append(line.split(":")[1]).append("\t");
						}
						return true;
					}

					@Override
					public String getResult() {
						return sb.toString();
					}

				});
		System.out.println(result);
	}
}
