package edu.cwru.cbc.ASM.commons.io;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.FASTQSequence;
import edu.cwru.cbc.ASM.commons.sequence.IUPACCode;

import java.io.IOException;
import java.util.LinkedHashSet;

/**
 * Created by lancelothk on 6/10/15.
 * LineProcessor for reading FASTQ format file.
 */
public class FASTQLineProcessor implements LineProcessor<LinkedHashSet<FASTQSequence>> {
	/**
	 * quality score range get from:
	 * Cock, Peter JA, et al. "The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants." Nucleic acids research 38.6 (2010): 1767-1771.
	 **/
	public static final int MIN_QUALITY_SCORE = 33;
	public static final int MAX_QUALITY_SCORE = 126;
	int lineCount = 0;
	// used to check if id is unique and also keep order of sequence as reading order.
	private LinkedHashSet<FASTQSequence> resultSet = new LinkedHashSet<>();
	private String id;
	private String sequence;

	@Override
	public boolean processLine(String line) throws IOException {
		switch (lineCount++) {
			case 0:
				processIdLine(line);
				break;
			case 1:
				processSequenceLine(line);
				break;
			case 2:
				processPlusLine(line);
				break;
			case 3:
				processQualityLine(line);
				break;
			default:
				throw new RuntimeException("incorrect line count!");
		}
		lineCount %= 4;
		return true;
	}

	private void processQualityLine(String line) {
		if (line.length() != sequence.length()) {
			throw new RuntimeException("length of quality String does not match length of sequence!");
		}
		/** There are different quality score systems. But all start from ASCII {@link #MIN_QUALITY_SCORE} and end before {@link #MAX_QUALITY_SCORE} **/
		for (char c : line.toCharArray()) {
			if (c < MIN_QUALITY_SCORE || c > MAX_QUALITY_SCORE) {
				throw new RuntimeException("invalid quality score character:" + c);
			}
		}
		if (!resultSet.add(new FASTQSequence(id, sequence, line))) {
			throw new RuntimeException("duplicate id detected:" + id);
		}
	}

	private void processPlusLine(String line) {
		if (line.charAt(0) != '+') {
			throw new RuntimeException("unexpected line!:" + line);
		}
		if (!line.substring(1, line.length()).equals(id)) {
			throw new RuntimeException(
					String.format("ID in plus line is not same to ID in ID line: plus:%s\tid:%s", line, id));
		}
	}

	private void processSequenceLine(String line) {
		if (!IUPACCode.validateNucleotideCode(line)) {
			throw new RuntimeException("invalid character in sequence!\t" + line);
		}
		sequence = line;
	}

	private void processIdLine(String line) {
		if (line.charAt(0) != '@') {
			throw new RuntimeException("unexpected line!:" + line);
		}
		if (line.length() == 1) {
			throw new RuntimeException("empty ID line!");
		}
		id = line.substring(1, line.length());
	}

	/**
	 * Result in LinkedHashSet keeps the sequences in the order of reading.
	 */
	@Override
	public LinkedHashSet<FASTQSequence> getResult() {
		switch (lineCount) {
			case 1:
				throw new RuntimeException("missing sequence line in the end of file!");
			case 2:
				throw new RuntimeException("missing plus line in the end of file!");
			case 3:
				throw new RuntimeException("missing quality line in the end of file!");
		}
		return resultSet;
	}
}
