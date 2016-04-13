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
	private RefChr refChr;
	private int rawIntervalCount;

	public CPMR(List<RefCpG> refCpGList, RefChr refChr, int min_cpg_coverage, int min_interval_cpg, int min_interval_reads, double partial_methyl_threshold) throws
			IOException {
		this.refChr = refChr;
		this.min_cpg_coverage = min_cpg_coverage;
		this.min_interval_cpg = min_interval_cpg;
		this.min_interval_reads = min_interval_reads;
		this.partial_methyl_threshold = partial_methyl_threshold;
		this.refCpGList = refCpGList;
	}

	public List<ImmutableGenomicInterval> getGenomicIntervals() {
		List<ImmutableGenomicInterval> immutableGenomicIntervalList = new ArrayList<>();
		// split regions by continuous CpG coverage
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
				if (resultRefCpGList.size() > 0) {
					filterInterval(resultRefCpGList, immutableGenomicIntervalList);

				}
				resultRefCpGList = new ArrayList<>();
			}
		}
		if (resultRefCpGList.size() > 0) {
			filterInterval(resultRefCpGList, immutableGenomicIntervalList);
		}
		return immutableGenomicIntervalList;
	}

	private void filterInterval(List<RefCpG> interval, List<ImmutableGenomicInterval> immutableGenomicIntervalList) {
		rawIntervalCount++;
		if (interval.size() >= min_interval_cpg) {
			// since reads are collected from refCpG, all reads contains at least one CpG.
			List<MappedRead> mappedReadList = getDistinctReadsFromRefCpGs(interval)
					.sorted(MappedRead::compareTo).collect(Collectors.toList());
			if (mappedReadList.size() >= min_interval_reads) {
				int startCpGPos = interval.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = interval.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				@SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
						int startPos = startCpGPos;
				int endPos = endCpGPos + 1;// +1 to include whole CpG in plus strand
				immutableGenomicIntervalList.add(
						new ImmutableGenomicInterval(refChr.getChr(), refChr.getRefString(), startPos, endPos, interval,
								mappedReadList));
			}
		}
	}

	private Stream<MappedRead> getDistinctReadsFromRefCpGs(List<RefCpG> list) {
		return list.stream().map(RefCpG::getCpGList).flatMap(Collection::stream).map(CpG::getMappedRead).distinct();
	}

	public int getRawIntervalCount() {
		return rawIntervalCount;
	}
}
