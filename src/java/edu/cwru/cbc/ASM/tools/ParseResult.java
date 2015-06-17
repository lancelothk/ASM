package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;

import java.io.File;
import java.io.IOException;

/**
 * Created by lancelothk on 6/16/15.
 * Used for parsing experiment intersection result.
 */
public class ParseResult {
	public static void main(String[] args) throws IOException {
		String result = Files.readLines(new File("/home/lancelothk/jun14_filterReads_intersection"), Charsets.UTF_8,
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
						} else if (line.startsWith("#covered regions in A")) {
							sb.append(line.split(":")[1].split("\\(")[0]).append("\n");
						} else if (line.startsWith("total length of intersection")) {
							sb.append(line.split(":")[1]).append("\t");
						}
						return true;
					}

					@Override
					public String getResult() {
						return sb.toString();
					}

				});
		Files.write(result, new File("/home/lancelothk/filter"), Charsets.UTF_8);
	}
}
