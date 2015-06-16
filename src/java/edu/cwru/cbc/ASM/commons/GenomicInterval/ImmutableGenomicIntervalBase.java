package edu.cwru.cbc.ASM.commons.GenomicInterval;

/**
 * Created by lancelothk on 6/8/15.
 * Immutable base class for genomic interval
 */
public class ImmutableGenomicIntervalBase {
	protected final String chr;
	protected final int start;
	protected final int end;

	public ImmutableGenomicIntervalBase(String chr, int start, int end) {
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
}
