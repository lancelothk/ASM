package edu.cwru.cbc.ASM.commons.genomicInterval;

/**
 * Created by kehu on 3/27/15.
 * Base class of Genomic interval.
 */
public abstract class GenomicIntervalBase implements Comparable<GenomicIntervalBase> {
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

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	@Override
	public int compareTo(GenomicIntervalBase o) {
		int chrCompare = this.chr.compareTo(o.chr);
		if (chrCompare != 0) {
			return chrCompare;
		}
		int startDiff = this.getStart() - o.getStart();
		int endDiff = this.getEnd() - o.getEnd();
		return startDiff == 0 ? endDiff : startDiff;
	}
}
