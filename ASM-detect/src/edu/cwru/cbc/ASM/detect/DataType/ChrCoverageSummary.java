package edu.cwru.cbc.ASM.detect.DataType;

import java.io.IOException;
import java.util.ArrayList;
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


	public ChrCoverageSummary(int chrSize) throws IOException {
		// start index is 1. Index range is 1 - chrSize
		this.chrBitSet = new BitSet(chrSize + 1);
		this.chrBitCoverage = new int[chrSize + 1];
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

	public List<GenomicInterval> generateIntervals() {
		List<GenomicInterval> genomicIntervals = new ArrayList<>();
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
			genomicIntervals.add(
					new GenomicInterval(setIndex, clearIndex - 1, clearIndex - setIndex, maxCount, avgCount));
		}
		return genomicIntervals;
	}
}
