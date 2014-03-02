package ASM.DataType;

import com.google.common.io.CharStreams;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;
import java.util.List;

/**
 * Created by lancelothk on 2/23/14.
 * Data structure for storing coverage status and count for each bit in chromosome
 */
// TODO add IntervalChekingLineProcessor into this class as an inner implementation of interface. If processor is used multiple places, it should keep separated.
public class ChrCoverageSummary {
	private BitSet chrBitSet;
	private int[] chrBitCoverage;
	private String reference;

	public ChrCoverageSummary(int chrSize, String referenceFileName) throws IOException {
		// start index is 1. Index range is 1 - chrSize
		this.chrBitSet = new BitSet(chrSize + 1);
		this.chrBitCoverage = new int[chrSize + 1];
		this.reference = readReference(referenceFileName);
	}

	public void addCoverage(int start, int end) {
		// BitSet.set exclude toIndex, so plus 1
		this.chrBitSet.set(start, end + 1);
		this.addBitCoverage(start, end);

	}

	private void addBitCoverage(int start, int end) {
		for (int i = start; i <= end; i++) {
			this.chrBitCoverage[i]++;
		}
	}

	public void writeIntervalCoverageSummary(String chr, BufferedWriter bufferedWriter) throws IOException {

		bufferedWriter.write("chr\tstart\tend\tlength\treadCount\tavgCount\tCpGCount\n");
		// bit of clearIndex is exclusive to interval
		int clearIndex = 0, setIndex = 0;
		while (clearIndex < chrBitSet.size() && setIndex < chrBitSet.size()) {
			setIndex = chrBitSet.nextSetBit(clearIndex);
			if (setIndex < 0) {
				break;
			}
			clearIndex = chrBitSet.nextClearBit(setIndex);
			int maxCount = 0;
			double avgCount = 0;
			for (int i = setIndex; i < clearIndex; i++) {
				if (chrBitCoverage[i] > maxCount) {
					maxCount = chrBitCoverage[i];
				}
				avgCount += chrBitCoverage[i];
			}
			avgCount = avgCount / (clearIndex - setIndex);
			int start = setIndex, end = clearIndex - 1;
			bufferedWriter.write(String.format("%s\t%d\t%d\t%d\t%d\t%f\t%d\n", chr, start, end, clearIndex - setIndex, maxCount, avgCount, countCpG(reference, start, end)));
		}
	}

	private String readReference(String fileName) throws IOException {
		List<String> lines = CharStreams.readLines(new BufferedReader(new FileReader(fileName)));
		StringBuilder referenceBuilder = new StringBuilder();
		lines.remove(0); // remove first line, which is chromosome name
		for (String line : lines) {
			referenceBuilder.append(line);
		}
		return referenceBuilder.toString();
	}

	private int countCpG(String reference, int start, int end) {
		// start and end is 1-based.
		int count = 0;
		for (int i = start - 1; i <= end - 1; i++) {
			if (i + 1 < reference.length() && (reference.charAt(i) == 'c' || reference.charAt(i) == 'C') && (reference.charAt(i + 1) == 'g' || reference.charAt(i + 1) == 'G')) { // 'C' is last character in reference
				count++;
			}
		}
		return count;
	}
}
