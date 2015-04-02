package edu.cwru.cbc.ASM.commons;

/**
 * Created by kehu on 3/27/15.
 * Base class of Genomic interval.
 */
public class GenomicIntervalBase {
	protected String chr;
	protected int start;
	protected int end;

	public GenomicIntervalBase(String chr, int start, int end) {
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
