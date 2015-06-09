package edu.cwru.cbc.ASM.commons.Read;

import edu.cwru.cbc.ASM.commons.CpG.RefCpG;

import java.io.IOException;
import java.util.BitSet;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * LineProcessor for process mapped read and print summary.
 */
public class MappedReadLineProcessorWithSummary extends MappedReadLineProcessorBase<String> {
	private long totalLength, qualifiedTotalLength;
	private long count;
	private int maxLength = Integer.MIN_VALUE;
	private int minLength = Integer.MAX_VALUE;
	private int qualifiedMaxLength = Integer.MIN_VALUE;
	private int qualifiedMinLength = Integer.MAX_VALUE;
	private long totalCpGCount;
	private long qualifiedTotalCpGCount;
	private int maxCpGCount = Integer.MIN_VALUE;
	private int minCpGCount = Integer.MAX_VALUE;
	private int qualifiedMaxCpGCount = Integer.MIN_VALUE;
	private int qualifiedMinCpGCount = Integer.MAX_VALUE;
	private BitSet chrBitSet;
	private BitSet informativeChrBitSet;
	private int min_read_cpg;

	public MappedReadLineProcessorWithSummary(List<RefCpG> refCpGList, int min_read_cpg, int refLength) throws
			IOException {
		super(refCpGList);
		this.min_read_cpg = min_read_cpg;
		this.chrBitSet = new BitSet(refLength);
		this.informativeChrBitSet = new BitSet(refLength);
	}

	@Override
	protected void updateRefCpG(MappedRead mappedRead) {
		int length = mappedRead.getSequence().length();
		int cpgCount = mappedRead.countCpG(refMap);
		if (cpgCount >= this.min_read_cpg) {
			super.updateRefCpG(mappedRead);
			qualifiedTotalLength += length;
			qualifiedMaxLength = length > qualifiedMaxLength ? length : qualifiedMaxLength;
			qualifiedMinLength = length < qualifiedMinLength ? length : qualifiedMinLength;
			qualifiedTotalCpGCount += cpgCount;
			qualifiedMaxCpGCount = cpgCount > qualifiedMaxCpGCount ? cpgCount : qualifiedMaxCpGCount;
			qualifiedMinCpGCount = cpgCount < qualifiedMinCpGCount ? cpgCount : qualifiedMinCpGCount;
			informativeChrBitSet.set(mappedRead.getStart(), mappedRead.getEnd(), true);
		}

		count++;
		totalLength += length;
		maxLength = length > maxLength ? length : maxLength;
		minLength = length < minLength ? length : minLength;
		totalCpGCount += cpgCount;
		maxCpGCount = cpgCount > maxCpGCount ? cpgCount : maxCpGCount;
		minCpGCount = cpgCount < minCpGCount ? cpgCount : minCpGCount;
		chrBitSet.set(mappedRead.getStart(), mappedRead.getEnd(), true);
	}

	public String getResult() {
		return getSummaryString();
	}

	private String getSummaryString() {
		long totalCpGCountInRefMap = 0;
		int maxCpGCoverage = Integer.MIN_VALUE;
		int minCpGCoverage = Integer.MAX_VALUE;
		for (RefCpG refCpG : refMap.values()) {
			int coverage = refCpG.getCpGCoverage();
			totalCpGCountInRefMap += coverage;
			maxCpGCoverage = coverage > maxCpGCoverage ? coverage : maxCpGCoverage;
			minCpGCoverage = coverage < minCpGCoverage ? coverage : minCpGCoverage;
		}
		StringBuilder sb = new StringBuilder();
		sb.append("Reference and Mapped Reads summary:\n");
		sb.append(String.format("Reference:\nrefLength:%d\trefCpgSiteNumber:%d\n", chrBitSet.size(), refMap.size()));
		sb.append(String.format(
				"\nRaw reads:\ntotalReadCount:%d\ttotalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n",
				count, totalLength, totalLength / (double) count, maxLength, minLength));
		sb.append(String.format("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n", totalCpGCount,
				totalCpGCount / (double) count, maxCpGCount, minCpGCount));
		sb.append(String.format("coveredLength:%d\toverallCoverage:%f\tactualCoverage:%f\n", chrBitSet.cardinality(),
				totalLength / (double) chrBitSet.size(), totalLength / (double) chrBitSet.cardinality()));

		sb.append(String.format(
				"\nInformative reads:\nreadCount:%d\ttotalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n",
				mappedReadList.size(), qualifiedTotalLength, qualifiedTotalLength / (double) mappedReadList.size(),
				qualifiedMaxLength, qualifiedMinLength));
		sb.append(String.format("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n",
				qualifiedTotalCpGCount, qualifiedTotalCpGCount / (double) mappedReadList.size(), qualifiedMaxCpGCount,
				qualifiedMinCpGCount));
		sb.append(String.format("coveredLength:%d\toverallCoverage:%f\tactualCoverage:%f\n",
				informativeChrBitSet.cardinality(),
				qualifiedTotalLength / (double) informativeChrBitSet.cardinality(),
				qualifiedTotalLength / (double) informativeChrBitSet.size()));
		sb.append(String.format(
				"(informative reads only)totalCpGCountInRefMap:%d\tavg cpgSite coverage:%f\tmaxCpGCoverage:%d\tminCpGCoverage:%d\n",
				totalCpGCountInRefMap, totalCpGCountInRefMap / (double) refMap.size(), maxCpGCoverage, minCpGCoverage));
		return sb.toString();
	}
}
