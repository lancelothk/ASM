package edu.cwru.cbc.ASM.commons.Read;

import edu.cwru.cbc.ASM.commons.CpG.RefCpG;

import java.io.IOException;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * LineProcessor for process mapped read and print summary.
 */
public class MappedReadLineProcessorWithSummary extends MappedReadLineProcessorBase<String> {
	private int refLength;
	private int min_read_cpg;
	private MappedReadSummary rawInputSummary;
	private MappedReadSummary filteredSummary;

	public MappedReadLineProcessorWithSummary(List<RefCpG> refCpGList, int min_read_cpg, int refLength) throws
			IOException {
		super(refCpGList);
		this.refLength = refLength;
		this.min_read_cpg = min_read_cpg;
		this.rawInputSummary = new MappedReadSummary(refLength);
		this.filteredSummary = new MappedReadSummary(refLength);
	}

	@Override
	protected void updateRefCpG(MappedRead mappedRead) {
		int cpgCount = mappedRead.countCpG(refMap);
		if (cpgCount >= this.min_read_cpg) {
			super.updateRefCpG(mappedRead);
			filteredSummary.addMappedRead(mappedRead, cpgCount);
		}
		rawInputSummary.addMappedRead(mappedRead, cpgCount);
	}

	public String getResult() {
		return getSummaryString();
	}

	private String getSummaryString() {
		return String.format("Reference:\nrefLength:%d\trefCpgSiteNumber:%d\n", refLength,
				refMap.size()) + rawInputSummary.getSummaryString(
				"\nSummary of raw input:\n") + filteredSummary.getSummaryString(
				"\nSummary of filtered reads:\n") + MappedReadSummary.getCpGCoverageSummary(
				refMap.values());
	}
}
