package edu.cwru.cbc.ASM.conversion.Conversion;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.apache.commons.cli.ParseException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by kehu on 2/6/15.
 * Convert Lister 2009 mapped reads to amrfinder required mr format.
 */
public class MappedReadToMR {

	public static void main(String[] args) throws IOException, ParseException {
		String input = "/home/lancelothk/experiments/ASM/data/i90_r1_chr20";
		String output = "/home/lancelothk/experiments/ASM/simulation/i90_r1_chr20.mr";
		conversion(input, output);
	}

	public static void conversion(String inputFileName, String outputFileName) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));

		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor() {

			@Override
			public boolean processLine(String line) throws IOException {
				try {
					if (line.startsWith("chr") || line.startsWith("ref")) {
						return true;
					} else if (line.equals("")) {
						return false;
					} else {
						String[] items = line.split("\t");
						if (items.length != 6) {
							throw new RuntimeException("invalid mapped read format:" + line);
						}
						if (items[1].length() != 1) {
							throw new RuntimeException("invalid strand!");
						}

						MappedRead mappedRead = new MappedRead(items[0], items[1].charAt(0), Integer.parseInt(items[2]),
								items[4], items[5]);
						writer.write(mappedRead.toMRFormatString(0, '~') + "\n");
						return true;
					}
				} catch (Exception e) {
					throw new RuntimeException("Problem line: " + line + "\n", e);
				}
			}

			@Override
			public Object getResult() {
				return null;
			}
		});
		writer.close();
	}

}
