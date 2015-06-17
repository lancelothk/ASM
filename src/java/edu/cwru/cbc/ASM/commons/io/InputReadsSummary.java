package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.util.BitSet;
import java.util.Collection;

/**
 * Created by lancelothk on 6/9/15.
 * Store and construct summary info about mapped reads.
 */
//TODO add unit test
public class InputReadsSummary {
	private long totalLength;
	private long count;
	private int maxLength = Integer.MIN_VALUE;
	private int minLength = Integer.MAX_VALUE;
	private long totalCpGCount;
	private int maxCpGCount = Integer.MIN_VALUE;
	private int minCpGCount = Integer.MAX_VALUE;
	private BitSet chrBitSet;
	private StringBuilder sb;
	private int refLength;

	public InputReadsSummary(int refLength) {
		this.refLength = refLength;
		this.sb = new StringBuilder();
		this.chrBitSet = new BitSet(refLength);
	}

	public static String getCpGCoverageSummary(Collection<RefCpG> refCpGList) {
		long totalCpGCountInRefMap = 0;
		int maxCpGCoverage = Integer.MIN_VALUE;
		int minCpGCoverage = Integer.MAX_VALUE;
		for (RefCpG refCpG : refCpGList) {
			int coverage = refCpG.getCpGCoverage();
			totalCpGCountInRefMap += coverage;
			maxCpGCoverage = coverage > maxCpGCoverage ? coverage : maxCpGCoverage;
			minCpGCoverage = coverage < minCpGCoverage ? coverage : minCpGCoverage;
		}
		return String.format(
				"#refCpGSite:%d\ttotalCpGCountInRefMap:%d\tavg cpgSite coverage:%f\tmaxCpGCoverage:%d\tminCpGCoverage:%d\n",
				refCpGList.size(), totalCpGCountInRefMap, totalCpGCountInRefMap / (double) refCpGList.size(),
				maxCpGCoverage, minCpGCoverage);
	}

	/**
	 * update summary with give read.
	 */
	public void addMappedRead(MappedRead mappedRead) {
		count++;
		int cpgCount = mappedRead.getCpgList().size();
		int length = mappedRead.getSequence().length();
		totalLength += length;
		maxLength = length > maxLength ? length : maxLength;
		minLength = length < minLength ? length : minLength;
		totalCpGCount += cpgCount;
		maxCpGCount = cpgCount > maxCpGCount ? cpgCount : maxCpGCount;
		minCpGCount = cpgCount < minCpGCount ? cpgCount : minCpGCount;
		chrBitSet.set(mappedRead.getStart(), mappedRead.getEnd(), true);
	}

	public String getSummaryString(String header) {
		sb.append(header);
		sb.append(String.format(
				"totalReadCount:%d\ttotalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n",
				count, totalLength, totalLength / (double) count, maxLength, minLength));
		sb.append(String.format("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n", totalCpGCount,
				totalCpGCount / (double) count, maxCpGCount, minCpGCount));
		sb.append(String.format("coveredLength:%d\toverallCoverage:%f\tactualCoverage:%f\n", chrBitSet.cardinality(),
				totalLength / (double) refLength, totalLength / (double) chrBitSet.cardinality()));
		return sb.toString();
	}
}
