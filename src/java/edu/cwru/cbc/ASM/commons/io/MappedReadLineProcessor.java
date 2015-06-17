package edu.cwru.cbc.ASM.commons.io;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by lancelothk on 6/8/15.
 * abstract base class for MappedReadLineProcessor.
 * Shared by both MappedReadLineProcessor and MappedReadLineProcessorWithSummary
 * Mapped reads start pos is 0-based, end pos is 0-based.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
	protected LinkedHashMap<String, MappedRead> mappedReadLinkedHashMap = new LinkedHashMap<>();

	// TODO check duplicates, keep order of reads as in original file.
	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr") || line.startsWith("ref") || line.startsWith("assembly")) {
			return true;
		} else {
			MappedRead mr = processRead(line);
			if (mappedReadLinkedHashMap.containsKey(mr.getId())) {
				throw new RuntimeException("found duplicate mapped read! in line:\t" + line);
			} else {
				mappedReadLinkedHashMap.put(mr.getId(), mr);
			}
			return true;
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

	private MappedRead processRead(String line) {
		String[] items = line.split("\t");
		if (!items[1].equals("+") && !items[1].equals("-")) {
			throw new RuntimeException("invalid strand! in line:\t" + line);
		}
		int start = Integer.parseInt(items[2]);// mapped read is 0 based start.

		if (Pattern.compile("[^ACGTN\\.]").matcher(items[4]).find()) {
			throw new RuntimeException("invalid character in sequence! only ACGTN and '.' are allowed!:\t" + line);
		}
		MappedRead mappedRead;
		if (items.length == 6) {
			// for h1/i90 dataset
			mappedRead = new MappedRead(items[0], items[1].charAt(0), start, items[4], items[5]);
		} else {
			throw new RuntimeException("columns is not correct for mapped read format in line:\t" + line);
		}
		return mappedRead;
	}

}
