package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.genomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by lancelothk on 6/8/15.
 * CPMR stands for Consecutively covered partial methylation region.
 */
public class CPMR {
	private int min_cpg_coverage;
	private int min_interval_cpg;
	private int min_interval_reads;
	private double partial_methyl_threshold;
	private List<RefCpG> refCpGList;
	private int rawIntervalCount;

	public CPMR(List<RefCpG> refCpGList, int min_cpg_coverage, int min_interval_cpg, int min_interval_reads, double partial_methyl_threshold) throws
			IOException {
		this.min_cpg_coverage = min_cpg_coverage;
		this.min_interval_cpg = min_interval_cpg;
		this.min_interval_reads = min_interval_reads;
		this.partial_methyl_threshold = partial_methyl_threshold;
		this.refCpGList = refCpGList;
		refCpGList.sort((RefCpG c1, RefCpG c2) -> c1.getPos() - c2.getPos());
		// assign cpg id after sorted
		for (int i = 0; i < refCpGList.size(); i++) {
			refCpGList.get(i).assignIndex(i);
		}
	}

	private List<List<RefCpG>> getRefCpGIntervals() {
		// split regions by continuous CpG coverage
		List<List<RefCpG>> cpgSiteIntervalList = new ArrayList<>();
		boolean cont = false;
		List<RefCpG> resultRefCpGList = new ArrayList<>();
		for (int i = 0; i < refCpGList.size() - 1; i++) {
			RefCpG curr = refCpGList.get(i);
			RefCpG next = refCpGList.get(i + 1);
			if (curr.hasCommonRead(next, min_cpg_coverage) && curr.hasPartialMethyl(
					partial_methyl_threshold) && next.hasPartialMethyl(partial_methyl_threshold)) {
				if (!cont) {
					resultRefCpGList.add(curr);
					cont = true;
				}
				resultRefCpGList.add(next);
			} else {
				cont = false;
				cpgSiteIntervalList.add(resultRefCpGList);
				resultRefCpGList = new ArrayList<>();
			}
		}
		if (resultRefCpGList.size() > 0) {
			cpgSiteIntervalList.add(resultRefCpGList);
		}
		rawIntervalCount = cpgSiteIntervalList.size();
		return cpgSiteIntervalList;
	}

	public List<ImmutableGenomicInterval> getGenomicIntervals(RefChr refChr) {
		List<ImmutableGenomicInterval> immutableGenomicIntervalList = new ArrayList<>();
		getRefCpGIntervals().stream().filter(list -> list.size() >= min_interval_cpg).forEach(list -> {
			// since reads are collected from refCpG, all reads contains at least one CpG.
			List<MappedRead> mappedReadList = getDistinctReadsFromRefCpGs(list)
					.sorted(MappedRead::compareTo).collect(Collectors.toList());
			if (mappedReadList.size() >= min_interval_reads) {
				int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				@SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
						int startPos = startCpGPos;
				int endPos = endCpGPos + 1;// +1 to include whole CpG in plus strand
				immutableGenomicIntervalList.add(
						new ImmutableGenomicInterval(refChr.getChr(), refChr.getRefString(), startPos, endPos, list,
								mappedReadList));
				}
		});
		return immutableGenomicIntervalList;
	}

	private Stream<MappedRead> getDistinctReadsFromRefCpGs(List<RefCpG> list) {
		return list.stream().map(RefCpG::getCpGList).flatMap(Collection::stream).map(CpG::getMappedRead).distinct();
	}

	public int getRawIntervalCount() {
		return rawIntervalCount;
	}
}
