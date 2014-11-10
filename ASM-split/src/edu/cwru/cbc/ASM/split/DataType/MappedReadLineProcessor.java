package edu.cwru.cbc.ASM.split.DataType;

import com.google.common.io.LineProcessor;

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
	private Map<Integer, CpGSite> refMap;
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
	private int filteredReads = 0;

	public MappedReadLineProcessor(Map<Integer, CpGSite> refMap, int refLength) {
		this.refMap = refMap;
		this.mappedReadList = new ArrayList<>();
		this.chrBitSet = new BitSet(refLength);
	}

	@Override
	public boolean processLine(String s) throws IOException {
		String[] items = s.split("\t");
		if (items.length != 6){
			throw new RuntimeException("invalid mapped read format:" + s);
		}
		// 0-chr;1-strand;2-start;3-end;4-sequence;5-id
		MappedRead mappedRead = new MappedRead(items[0], items[1], items[5]);
		//to suit prev program
		mappedRead.setSequence(items[4]);
		mappedRead.setStart(Integer.parseInt(items[2]));
		mappedRead.setEnd(Integer.parseInt(items[3]));

		int start = Integer.parseInt(items[2]);// mapped read is 1bp right than UCSC ref
		int end = Integer.parseInt(items[3]);

		for (int i = start; i < end; i++) {
			if (refMap.containsKey(i)){
				CpG cpg = new CpG(mappedRead, refMap.get(i));
				cpg.setMethylated(extractMethylStatus(items[1], items[4].substring(i - start, i - start + 2)));
				mappedRead.addCpG(cpg);
				i++;
			}
		}

		int length = mappedRead.getSequence().length();
		int cpgCount = mappedRead.getCpGCount();

		// only use the read if contains at least two CpG sites
		if (cpgCount > 1) {
			mappedReadList.add(mappedRead);
			for (CpG cpg : mappedRead.getCpgList()) {
				refMap.get(cpg.getPos()).addCpG(cpg);
			}
			filteredReads++;
		}

		//stat variables
		count++;
		totalLength += length;
		maxLength = length > maxLength?length:maxLength;
		minLength = length < minLength?length:minLength;
		totalCpGCount += cpgCount;
		maxCpGCount = cpgCount > maxCpGCount?cpgCount:maxCpGCount;
		minCpGCount = cpgCount < minCpGCount?cpgCount:minCpGCount;
		chrBitSet.set(mappedRead.getStart()-1, mappedRead.getEnd() -1, true);
		countCoverCpG += (cpgCount > 0? 1:0);
		return true;
	}

	private boolean extractMethylStatus(String strand, String sequence) {
		if (sequence.length() == 1){
			return false;
		}
		else if (sequence.length() > 2){
			throw new RuntimeException("invalid input cpg sequence:\t" + sequence);
		}
		// only consider 'C' position in CpG. Ignore the char in 'G' position
		switch (strand) {
			case "+": {
				// CN is methylated, TN is non-methylated
				char cbp = sequence.charAt(0);
				char gbp = sequence.charAt(1);
				if (cbp == 'C' && gbp == 'G') {
					return true;
//				} else if (bp == 'T') {
//					return false;
				} else {
					// TODO should ignore non C/T case in CpG list
					return false;
				}
			}
			case "-": {
				// NC is methylated, NT is non-methylated
				char cbp = sequence.charAt(1);
				char gbp = sequence.charAt(0);
				if (cbp == 'C' && gbp == 'G') {  // for complementary bp, 'G'
					return true;
//				} else if (bp == 'T') { // for complementary bp, 'A'
//					return false;
				} else {
					// TODO should ignore non C/T case in CpG list
					return false;
				}
			}
			default:
				throw new RuntimeException("illegal strand!");
		}
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
		System.out.printf("totalCpGCount:%d\tavgCpgCount:%f\tmaxCpgCount:%d\tminCpgCount:%d\n",
						  totalCpGCount, totalCpGCount / (double) count, maxCpGCount, minCpGCount);
		System.out.printf("refLength:%d\trefCpgSiteNumber:%d\n", chrBitSet.size(), refMap.size());
		System.out.printf("readCoverage:%f\tcoveredLength:%d\n",
						  totalLength / (double) chrBitSet.size(),  chrBitSet.cardinality());
		System.out.println("FilteredReads:\t" + filteredReads);
		long totalCpGCountInRefMap = 0;
		int maxCpGCoverage = Integer.MIN_VALUE;
		int minCpGCoverage = Integer.MAX_VALUE;
		for (CpGSite cpGSite : refMap.values()) {
			int coverage = cpGSite.getCoverage();
			totalCpGCountInRefMap += coverage;
			maxCpGCoverage = coverage > maxCpGCoverage?coverage:maxCpGCoverage;
			minCpGCoverage = coverage < minCpGCoverage?coverage:minCpGCoverage;
		}
		System.out.printf("totalCpGCountInRefMap:%d\tavg cpgSite coverage:%f\tmaxCpGCoverage:%d\tminCpGCoverage:%d\n",
						  totalCpGCountInRefMap, totalCpGCountInRefMap / (double) refMap.size(), maxCpGCoverage, minCpGCoverage);
	}
}
