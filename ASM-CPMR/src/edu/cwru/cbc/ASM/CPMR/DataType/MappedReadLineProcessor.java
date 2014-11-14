package edu.cwru.cbc.ASM.CPMR.DataType;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MethylStatus;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

/**
 * Created by lancelothk on 5/27/14.
 * LineProcessor for process mapped read.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
	private Map<Integer, RefCpG> refMap;
	private List<MappedRead> mappedReadList;
	private long totalLength = 0;
	private long count = 0;
	private int maxLength = Integer.MIN_VALUE;
	private int minLength = Integer.MAX_VALUE;
	private long totalCpGCount = 0;
	private int maxCpGCount = Integer.MIN_VALUE;
	private int minCpGCount = Integer.MAX_VALUE;
	private BitSet chrBitSet;
	private int countCoverCpG = 0;

	public MappedReadLineProcessor(Map<Integer, RefCpG> refMap, int refLength) {
		this.refMap = refMap;
		this.mappedReadList = new ArrayList<>();
		this.chrBitSet = new BitSet(refLength);
	}

	@Override
	public boolean processLine(String s) throws IOException {
		String[] items = s.split("\t");
		if (items.length != 6) {
			throw new RuntimeException("invalid mapped read format:" + s);
		}
		// 0-chr;1-strand;2-start;3-end;4-sequence;5-id
		if (items[1].length() != 1) {
			throw new RuntimeException("invalid strand!");
		}
		int start = Integer.parseInt(items[2]);// mapped read is 1bp right than UCSC ref
		int end = Integer.parseInt(items[3]);
		MappedRead mappedRead = new MappedRead(items[0], items[1].charAt(0), start, end, items[4], items[5]);
		for (int i = start; i < end; i++) {
			if (refMap.containsKey(i)) {
				CpG cpg = new CpG(mappedRead, refMap.get(i));
				cpg.setMethylStatus(mappedRead.getMethylStatus(i));
				// ignore unknown methyl CpG
				if (cpg.getMethylStatus() != MethylStatus.N) {
					mappedRead.addCpG(cpg);
				}
				i++;
			}
		}

		int length = mappedRead.getSequence().length();
		int cpgCount = mappedRead.getCpGCount();

		mappedReadList.add(mappedRead);
		for (CpG cpg : mappedRead.getCpgList()) {
			refMap.get(cpg.getPos()).addCpG(cpg);
		}

		//stat variables
		count++;
		totalLength += length;
		maxLength = length > maxLength ? length : maxLength;
		minLength = length < minLength ? length : minLength;
		totalCpGCount += cpgCount;
		maxCpGCount = cpgCount > maxCpGCount ? cpgCount : maxCpGCount;
		minCpGCount = cpgCount < minCpGCount ? cpgCount : minCpGCount;
		chrBitSet.set(mappedRead.getStart() - 1, mappedRead.getEnd() - 1, true);
		countCoverCpG += (cpgCount > 0 ? 1 : 0);
		return true;
	}

	@Override
	public List<MappedRead> getResult() {
		printSummary();
		return this.mappedReadList;
	}

	private void printSummary() {
		System.out.println("Reference and Mapped Reads summary:");
		System.out.printf("totalReadCount:%d\tcountCoverAtLeastOneCpG:%d\n", count, countCoverCpG);
		System.out.printf("totalLength:%d\tavgLength:%f\tmaxLength:%d\tminLength:%d\n", totalLength,
						  totalLength / (double) count, maxLength, minLength);
		System.out.printf("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n", totalCpGCount,
						  totalCpGCount / (double) count, maxCpGCount, minCpGCount);
		System.out.printf("refLength:%d\trefCpgSiteNumber:%d\n", chrBitSet.size(), refMap.size());
		System.out.printf("readCoverage:%f\tcoveredLength:%d\n", totalLength / (double) chrBitSet.size(),
						  chrBitSet.cardinality());
		long totalCpGCountInRefMap = 0;
		int maxCpGCoverage = Integer.MIN_VALUE;
		int minCpGCoverage = Integer.MAX_VALUE;
		for (RefCpG refCpG : refMap.values()) {
			int coverage = refCpG.getCoverage();
			totalCpGCountInRefMap += coverage;
			maxCpGCoverage = coverage > maxCpGCoverage ? coverage : maxCpGCoverage;
			minCpGCoverage = coverage < minCpGCoverage ? coverage : minCpGCoverage;
		}
		System.out.printf("totalCpGCountInRefMap:%d\tavg cpgSite coverage:%f\tmaxCpGCoverage:%d\tminCpGCoverage:%d\n",
						  totalCpGCountInRefMap, totalCpGCountInRefMap / (double) refMap.size(), maxCpGCoverage,
						  minCpGCoverage);
	}
}
