package edu.cwru.cbc.ASM.commons.DataType;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.BitSet;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * LineProcessor for process mapped read and print summary.
 */
public class MappedReadLineProcessorWithSummary extends MappedReadLineProcessor {
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
	private BufferedWriter summaryWriter;

	public MappedReadLineProcessorWithSummary(List<RefCpG> refCpGList, int min_read_cpg, int refLength,
											  String reportFileName) throws IOException {
		super(refCpGList, min_read_cpg);
		this.chrBitSet = new BitSet(refLength);
		this.summaryWriter = new BufferedWriter(new FileWriter(reportFileName));
	}

	@Override
	public boolean processLine(String s) throws IOException {
		try {
			MappedRead mappedRead = processRead(s);
			int length = mappedRead.getSequence().length();
			int cpgCount = mappedRead.getCpgList().size();

			//stat variables
			if (cpgCount >= this.min_read_cpg) {
				qualifiedTotalLength += length;
				qualifiedMaxLength = length > qualifiedMaxLength ? length : qualifiedMaxLength;
				qualifiedMinLength = length < qualifiedMinLength ? length : qualifiedMinLength;
				qualifiedTotalCpGCount += cpgCount;
				qualifiedMaxCpGCount = cpgCount > qualifiedMaxCpGCount ? cpgCount : qualifiedMaxCpGCount;
				qualifiedMinCpGCount = cpgCount < qualifiedMinCpGCount ? cpgCount : qualifiedMinCpGCount;
			}
			count++;
			totalLength += length;
			maxLength = length > maxLength ? length : maxLength;
			minLength = length < minLength ? length : minLength;
			totalCpGCount += cpgCount;
			maxCpGCount = cpgCount > maxCpGCount ? cpgCount : maxCpGCount;
			minCpGCount = cpgCount < minCpGCount ? cpgCount : minCpGCount;
			chrBitSet.set(mappedRead.getStart() - 1, mappedRead.getEnd() - 1, true);
			return true;
		} catch (Exception e) {
			throw new RuntimeException("Problem line: " + s + "\n", e);
		}
	}

	@Override
	public List<MappedRead> getResult() {
		try {
			printSummary();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return this.mappedReadList;
	}

	private void printSummary() throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("Reference and Mapped Reads summary:\n");
		sb.append(String.format("totalReadCount:%d\tQualifiedReadCount:%d\n", count, mappedReadList.size()));
		sb.append(String.format("totalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n", totalLength,
				totalLength / (double) count, maxLength, minLength));
		sb.append(String.format("qualified: totalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n",
				qualifiedTotalLength, qualifiedTotalLength / (double) mappedReadList.size(), qualifiedMaxLength,
				qualifiedMinLength));
		sb.append(String.format("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n", totalCpGCount,
				totalCpGCount / (double) count, maxCpGCount, minCpGCount));
		sb.append(String.format("qualified: totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n",
				qualifiedTotalCpGCount, qualifiedTotalCpGCount / (double) mappedReadList.size(), qualifiedMaxCpGCount,
				qualifiedMinCpGCount));
		sb.append(String.format("refLength:%d\trefCpgSiteNumber:%d\n", chrBitSet.size(), refMap.size()));
		sb.append(String.format("readCoverage:%f\tcoveredLength:%d\n", totalLength / (double) chrBitSet.size(),
				chrBitSet.cardinality()));
		long totalCpGCountInRefMap = 0;
		int maxCpGCoverage = Integer.MIN_VALUE;
		int minCpGCoverage = Integer.MAX_VALUE;
		for (RefCpG refCpG : refMap.values()) {
			int coverage = refCpG.getCpGCoverage();
			totalCpGCountInRefMap += coverage;
			maxCpGCoverage = coverage > maxCpGCoverage ? coverage : maxCpGCoverage;
			minCpGCoverage = coverage < minCpGCoverage ? coverage : minCpGCoverage;
		}
		sb.append(String.format(
				"totalCpGCountInRefMap:%d\tavg cpgSite coverage:%f\tmaxCpGCoverage:%d\tminCpGCoverage:%d\n",
				totalCpGCountInRefMap, totalCpGCountInRefMap / (double) refMap.size(), maxCpGCoverage, minCpGCoverage));
		summaryWriter.write(sb.toString());
		summaryWriter.close();
	}
}
