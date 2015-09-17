package edu.cwru.cbc.ASM.commons.methylation;

/**
 * Created by kehu on 9/17/15.
 */
public class RefCpGStat {
	private final String chr;
	private final int pos;
	private final int coveredCount;
	private final double methylLevel;

	public RefCpGStat(String chr, int pos, int coveredCount, double methylLevel) {
		this.chr = chr;
		this.pos = pos;
		this.coveredCount = coveredCount;
		this.methylLevel = methylLevel;
	}

	public int getPos() {
		return pos;
	}

	public int getCoveredCount() {
		return coveredCount;
	}

	public double getMethylLevel() {
		return methylLevel;
	}
}
