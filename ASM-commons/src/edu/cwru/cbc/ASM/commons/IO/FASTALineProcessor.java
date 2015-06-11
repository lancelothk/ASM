package edu.cwru.cbc.ASM.commons.IO;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.Sequence.FASTASequence;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

/**
 * Created by lancelothk on 6/10/15.
 * LineProcessor for reading FASTA format file.
 */
public class FASTALineProcessor implements LineProcessor<LinkedHashMap<String, FASTASequence>> {
	// used to check if id is unique and also keep order of sequence as reading order.
	private LinkedHashMap<String, FASTASequence> resultMap = new LinkedHashMap<>();
	private int stringBuilderInitCapacity = -1; // -1 means use default setting.
	private String id = null;
	private StringBuilder sb;

	public FASTALineProcessor() {
		this.sb = new StringBuilder();
	}

	/**
	 * For reading ultra-long sequence. E.g. reference sequence.
	 */
	public FASTALineProcessor(int stringBuilderInitCapacity) {
		this.stringBuilderInitCapacity = stringBuilderInitCapacity;
		// constructor of StringBuilder will check if the capacity value is valid. No need to do it here.
		this.sb = new StringBuilder(stringBuilderInitCapacity);
	}

	/**
	 * 1. Support both single/multi fasta format.
	 * 2. When meet multiple sequence lines, combine them into one line.
	 * 3. Sequence should not contain character other than acgtnACGTN and '.'.
	 * Notice: this method shouldn't be used individually without calling {@link #getResult()} to finish the file
	 * processing.
	 */
	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith(">")) {
			if (line.length() == 1) {
				throw new RuntimeException("empty id!");
			}
			if (id == null) { // first sequence
				id = line.substring(1, line.length());
			} else {
				if (sb.length() == 0) {
					throw new RuntimeException("missing sequence line!");
				}
				if (resultMap.containsKey(id)) {
					throw new RuntimeException("duplicate id detected!");
				} else if (sb.length() != 0) {
					resultMap.put(id, new FASTASequence(id, sb.toString()));
					if (stringBuilderInitCapacity == -1) {
						sb = new StringBuilder();
					} else {
						sb = new StringBuilder(stringBuilderInitCapacity);
					}
				}
				id = line.substring(1, line.length());
			}
		} else {
			if (id == null) {
				throw new RuntimeException("missing '>' in the beginning of file!");
			}
			if (Pattern.compile("[^ACGTNacgtn\\.]").matcher(line).find()) {
				throw new RuntimeException("invalid character in sequence! only acgtnACGTN and '.' are allowed!");
			}
			sb.append(line);
		}
		return true;
	}

	@Override
	public LinkedHashMap<String, FASTASequence> getResult() {
		// put last FASTA sequence to map.
		if (id != null) {
			if (sb.length() == 0) {
				throw new RuntimeException("missing sequence line!");
			} else {
				if (resultMap.containsKey(id)) {
					throw new RuntimeException("duplicate id detected!");
				} else {
					resultMap.put(id, new FASTASequence(id, sb.toString()));
				}
			}
		}
		return resultMap;
	}
}
