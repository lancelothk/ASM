package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;

import java.util.BitSet;
import java.util.Collection;

/**
 * Created by lancelothk on 6/9/15.
 * Store and construct summary info about mapped reads.
 */
public class InputReadsSummary {
	private long totalReadLength;
	private long readCount;
	private int maxReadLength = Integer.MIN_VALUE;
	private int minReadLength = Integer.MAX_VALUE;
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
		readCount++;
		int cpgCount = mappedRead.getCpgList().size();
		int readLength = mappedRead.getSequence().length();
		totalReadLength += readLength;
		maxReadLength = readLength > maxReadLength ? readLength : maxReadLength;
		minReadLength = readLength < minReadLength ? readLength : minReadLength;
		totalCpGCount += cpgCount;
		maxCpGCount = cpgCount > maxCpGCount ? cpgCount : maxCpGCount;
		minCpGCount = cpgCount < minCpGCount ? cpgCount : minCpGCount;
		chrBitSet.set(mappedRead.getStart(), mappedRead.getEnd(), true);
	}

	public String getSummaryString(String header) {
		sb.append(header);
		sb.append(String.format(
				"totalReadCount:%d\ttotalReadLength:%d\tavgLength:%f\tmaxReadLength:%d\tminReadLength:%d\n",
				readCount, totalReadLength, totalReadLength / (double) readCount, maxReadLength, minReadLength));
		sb.append(String.format("referenceLength:%d\tcoveredLength:%d\toverallCoverage:%f\tactualCoverage:%f\n",
				refLength, chrBitSet.cardinality(), totalReadLength / (double) refLength,
				totalReadLength / (double) chrBitSet.cardinality()));
		sb.append(
				String.format("totalCpGCount:%d\tavgCpgCountPerRead:%f\tmaxCpgCountPerRead:%d\tminCpgCountPerRead:%d\n",
						totalCpGCount, totalCpGCount / (double) readCount, maxCpGCount, minCpGCount));
		return sb.toString();
	}

	public String getSummaryString(String header, HashIntObjMap<RefCpG> refMap) {
		long coveredCpGSites = refMap.values().stream().filter(cpg -> cpg.getCoveredCount() > 0).count();
		return getSummaryString(header) + (String.format(
				"ReferenceCpGSitesCount:%d\tCoveredCpGSites:%d\tCpGSiteCoverage:%f\tActualCpGSiteCoverage:%f\n",
				refMap.size(), coveredCpGSites, totalCpGCount / (double) refMap.size(),
				totalCpGCount / (double) coveredCpGSites));
	}
}
