package edu.cwru.cbc.ASM.tools.conversion;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by kehu on 2/6/15.
 * convert amrfinder MR format to Lister 2009 mapped reads format.
 */
public class MRToMappedRead {
	public static void main(String[] args) throws IOException {
		String input = "/home/kehu/experiments/ASM/amrfinder/simulation/i90_r1_chr20_sim.mr";
		String output = "/home/kehu/experiments/ASM/amrfinder/simulation/i90_r1_chr20_sim";
		conversion(input, output);
	}

	public static void conversion(String input, String output) throws IOException {
		BufferedWriter mappedReadWriter = new BufferedWriter(new FileWriter(output));

		Files.readLines(new File(input), Charsets.UTF_8, new LineProcessor() {

			@Override
			public boolean processLine(String line) throws IOException {
				String[] items = line.split("\t");
				if (items.length != 8) {
					throw new RuntimeException("invalid mapped read format:" + line);
				}
				if (items[5].length() != 1) {
					throw new RuntimeException("invalid strand!");
				}

				MappedRead mappedRead = new MappedRead(items[0], items[5].charAt(0), Integer.parseInt(items[1]),
						items[6], Integer.parseInt(items[3]));
				mappedReadWriter.write(mappedRead.toMappedReadString());
				return true;
			}

			@Override
			public Object getResult() {
				return null;
			}

		});


		mappedReadWriter.close();
	}
}
