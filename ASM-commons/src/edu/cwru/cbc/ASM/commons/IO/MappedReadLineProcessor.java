package edu.cwru.cbc.ASM.commons.IO;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by lancelothk on 6/8/15.
 * abstract base class for MappedReadLineProcessor.
 * Shared by both MappedReadLineProcessor and MappedReadLineProcessorWithSummary
 * Mapped reads start pos is 0-based, end pos is 0-based.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
	protected List<MappedRead> mappedReadList = new ArrayList<>();

	// TODO check duplicates, keep order of reads as in original file.
	@Override
	public boolean processLine(String line) throws IOException {
		try {
			if (line.startsWith("chr") || line.startsWith("ref") || line.startsWith("assembly")) {
				return true;
			} else {
				mappedReadList.add(processRead(line));
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

	protected MappedRead processRead(String line) {
		String[] items = line.split("\t");
		if (!items[1].equals("+") && !items[1].equals("-")) {
			throw new RuntimeException("invalid strand!");
		}
		int start = Integer.parseInt(items[2]);// mapped read is 0 based start.

		if (Pattern.compile("[^ACGTN\\.]").matcher(items[4]).find()) {
			throw new RuntimeException("invalid character in sequence! only ACGTN and '.' are allowed!:\t" + line);
		}
		MappedRead mappedRead;
		if (items.length == 6) {
			// for h1/i90 dataset
			mappedRead = new MappedRead(items[0], items[1].charAt(0), start, items[4], items[5]);
		} else if (items.length == 7) {
			// for h9 and other dataset
			mappedRead = new MappedRead(items[0], items[1].charAt(0), start, items[4], items[6]);
		} else {
			throw new RuntimeException("columns is not correct for mapped read format");
		}
		return mappedRead;
	}

}
