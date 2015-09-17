package edu.cwru.cbc.ASM.commons.genomicInterval;

/**
 * Created by kehu on 9/17/15.
 * BedInterval for MethylStat query
 */
public class MethylStatBedInterval extends GenomicIntervalBase {
	private int coveredCount = 0;
	private int totalCount = 0;
	private double sumCoverage = 0;
	private double sumMethylLevel = 0;
	private double weightedSumMethylLevel = 0;
	private String name;

	public MethylStatBedInterval(String chr, int start, int end, String name) {
		super(chr, start, end);
		this.name = name;
	}

	public void intersect(String chr, int pos, int coverage, double methylLevel) {
		if (this.chr.equals(chr) && pos >= this.start && pos <= this.end) {
			if (coverage > 0) {
				this.coveredCount++;
				this.sumCoverage += coverage;
				this.sumMethylLevel += methylLevel;
				this.weightedSumMethylLevel += coverage * methylLevel;
			}
			this.totalCount++;
		}
	}

	@Override
	public String toString() {
		return String.format("%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f", this.name, this.chr, this.start, this.end,
				totalCount, coveredCount, sumCoverage / coveredCount, sumCoverage / totalCount,
				sumMethylLevel / coveredCount, sumMethylLevel / totalCount, weightedSumMethylLevel / sumCoverage);
	}
}
