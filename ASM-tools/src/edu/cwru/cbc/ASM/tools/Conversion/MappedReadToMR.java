package edu.cwru.cbc.ASM.tools.Conversion;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 2/6/15.
 * Convert Lister 2009 mapped reads to amrfinder required mr format.
 */
public class MappedReadToMR {

	public static void main(String[] args) throws IOException {
		String input = "/home/kehu/experiments/ASM/data/i90_r1_chr20";
		String output = "/home/kehu/experiments/ASM/data/i90_r1_chr20.mr";

		List<MappedRead> mappedReadList = Files.readLines(new File(input), Charsets.UTF_8,
				new LineProcessor<List<MappedRead>>() {
					private List<MappedRead> mappedReadList = new ArrayList<>();

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

								MappedRead mappedRead = new MappedRead(items[0], items[1].charAt(0),
										Integer.parseInt(items[2]), Integer.parseInt(items[3]), items[4], items[5]);
								mappedReadList.add(mappedRead);
								return true;
							}
						} catch (Exception e) {
							throw new RuntimeException("Problem line: " + line + "\n", e);
						}
					}

					@Override
					public List<MappedRead> getResult() {
						return mappedReadList;
					}
				});

		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		for (MappedRead mappedRead : mappedReadList) {
			// '~' is highest quality score.
			writer.write(mappedRead.toMRFormatString(0, '~') + "\n");
		}
		writer.close();

	}

}
