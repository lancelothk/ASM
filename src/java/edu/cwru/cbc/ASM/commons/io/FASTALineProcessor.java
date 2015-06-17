package edu.cwru.cbc.ASM.commons.io;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.FASTASequence;

import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.regex.Pattern;

/**
 * Created by lancelothk on 6/10/15.
 * LineProcessor for reading FASTA format file.
 */
public class FASTALineProcessor implements LineProcessor<LinkedHashSet<FASTASequence>> {
	// used to check if id is unique and also keep order of sequence as reading order.
	private LinkedHashSet<FASTASequence> resultSet = new LinkedHashSet<>();
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
		if (line.charAt(0) == '>') {
			if (line.length() == 1) {
				throw new RuntimeException("empty id!");
			}
			if (id == null) { // first sequence
				id = line.substring(1, line.length());
			} else {
				processFASTASequence();
				sb = initializeStringBuilder();
				id = line.substring(1, line.length());
			}
		} else {
			if (id == null) {
				throw new RuntimeException("missing '>' in the beginning of file!");
			}
			if (Pattern.compile("[^ACGTNacgtn\\.]").matcher(line).find()) {
				throw new RuntimeException("invalid character in sequence! only acgtnACGTN and '.' are allowed!:\t" + line);
			}
			sb.append(line);
		}
		return true;
	}

	private StringBuilder initializeStringBuilder() {
		if (stringBuilderInitCapacity == -1) {
			return new StringBuilder();
		} else {
			return new StringBuilder(stringBuilderInitCapacity);
		}
	}

	/**
	 * Result in LinkedHashSet keeps the sequences in the order of reading.
	 */
	@Override
	public LinkedHashSet<FASTASequence> getResult() {
		// put last FASTA sequence to map.
		if (id != null) {
			processFASTASequence();
		}
		return resultSet;
	}

	private void processFASTASequence() {
		if (sb.length() == 0) {
			throw new RuntimeException("missing sequence line!");
		}
		if (!resultSet.add(new FASTASequence(id, sb.toString()))) {
			throw new RuntimeException("duplicate id:" + id);
		}
	}
}
