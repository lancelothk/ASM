package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.GenomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.Methylation.CpG;
import edu.cwru.cbc.ASM.commons.Methylation.RefChr;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by lancelothk on 6/8/15.
 * CPMR stands for Consecutive covered partial methylation region.
 */
public class CPMR {
	private int min_cpg_coverage;
	private int min_read_cpg;
	private int min_interval_cpg;
	private int min_interval_reads;
	private double partial_methyl_threshold;
	private List<RefCpG> refCpGList;

	public CPMR(List<RefCpG> refCpGList, int min_cpg_coverage, int min_read_cpg,
	            int min_interval_cpg, int min_interval_reads, double partial_methyl_threshold) throws IOException {
		this.min_cpg_coverage = min_cpg_coverage;
		this.min_read_cpg = min_read_cpg;
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

	public List<List<RefCpG>> getIntervals() {
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
		return cpgSiteIntervalList;
	}

	public List<ImmutableGenomicInterval> getGenomicIntervals(RefChr refChr, List<List<RefCpG>> cpgSiteIntervalList) {
		List<ImmutableGenomicInterval> immutableGenomicIntervalList = new ArrayList<>();
		cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = getReadsFromRefCpGs(min_cpg_coverage, list);
			// only pass high quality result for next step.
			if (mappedReadSet.size() >= min_interval_reads && list.size() >= min_interval_cpg) {
				List<MappedRead> mappedReadList = mappedReadSet.stream().filter(
						mr -> countOverlapCpG(mr, list) >= min_read_cpg)
						.sorted((m1, m2) -> m1.getId().compareTo(m2.getId()))
						.sorted((m1, m2) -> m1.getStart() - m2.getStart())
						.collect(Collectors.toList());
				int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				@SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
						int startPos = startCpGPos;
				int endPos = endCpGPos + 1;
				immutableGenomicIntervalList.add(
						new ImmutableGenomicInterval(refChr.getChr(), refChr.getRefString(), startPos, endPos, list,
								mappedReadList));
			}
		});
		return immutableGenomicIntervalList;
	}

	private int countOverlapCpG(MappedRead read, List<RefCpG> refCpGList) {
		Map<Integer, RefCpG> refCpGMap = refCpGList.stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		int count = 0;
		for (CpG cpG : read.getCpgList()) {
			if (refCpGMap.containsKey(cpG.getPos())) {
				count++;
			}
		}
		return count;
	}

	private Set<MappedRead> getReadsFromRefCpGs(int min_cpg_coverage, List<RefCpG> list) {
		Set<MappedRead> mappedReadSet = new HashSet<>();
		list.forEach((refCpG) -> {
			refCpG.getCpGList().forEach((CpG cpg) -> mappedReadSet.add(cpg.getMappedRead()));
			assert refCpG.getCpGList().size() >= min_cpg_coverage; // make sure coverage is correct
		});
		return mappedReadSet;
	}
}
