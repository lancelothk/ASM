package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.Constant;
import edu.cwru.cbc.ASM.commons.CpG.CpG;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by lancelothk on 6/8/15.
 * CPMR stands for Consecutive covered partial methylation region.
 */
public class CPMR {
	public static final int INIT_POS = 0;
	private static int outputIntervalCount = 0;
	private int min_cpg_coverage;
	private int min_read_cpg;
	private int min_interval_cpg;
	private int min_interval_reads;
	private List<RefCpG> refCpGList;
	private RefChr refChr;

	public CPMR(RefChr refChr, List<RefCpG> refCpGList, int min_cpg_coverage, int min_read_cpg,
	            int min_interval_cpg, int min_interval_reads) throws IOException {
		this.min_cpg_coverage = min_cpg_coverage;
		this.min_read_cpg = min_read_cpg;
		this.min_interval_cpg = min_interval_cpg;
		this.min_interval_reads = min_interval_reads;
		this.refChr = refChr;
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
			if (curr.hasCommonRead(next, min_cpg_coverage) && curr.hasPartialMethyl() && next.hasPartialMethyl()) {
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

	public void writeIntervals(String intervalFolderName, String summaryFileName,
	                           List<List<RefCpG>> cpgSiteIntervalList) throws IOException {
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
		intervalSummaryWriter.write("chr\tlength\treadCount\tCpGCount\tstartCpG\tendCpG\n");
		cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = getReadsFromCpGs(min_cpg_coverage, list);
			// only pass high quality result for next step.
			if (mappedReadSet.size() >= min_interval_reads && list.size() >= min_interval_cpg) {
				List<MappedRead> mappedReadList = mappedReadSet.stream().filter(
						mr -> countOverlapCpG(mr, list) >= min_read_cpg)
						.sorted((m1, m2) -> m1.getId().compareTo(m2.getId()))
						.sorted((m1, m2) -> m1.getStart() - m2.getStart())
						.collect(Collectors.toList());
				outputIntervalCount++;
				int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				@SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
						int startPos = startCpGPos;
				int endPos = endCpGPos + 1;
				try {
					intervalSummaryWriter.write(
							String.format("%s\t%d\t%d\t%d\t%d\t%d\n", refChr.getChr(), endPos - startPos + 1,
									mappedReadList.size(), list.size(), startCpGPos, endCpGPos));
					writeMappedReadInInterval(intervalFolderName, refChr, startPos, endPos, mappedReadList);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		});
		intervalSummaryWriter.close();
	}

	public void writeReport(String reportFileName, String reportString, List<List<RefCpG>> cpgSiteIntervalList) throws
			IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(reportFileName, true));
		writer.write(reportString);
		writer.write("Raw Interval count:\t" + cpgSiteIntervalList.size() + "\n");
		writer.write("Output Interval count:\t" + outputIntervalCount + "\n");
		writer.close();
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

	private Set<MappedRead> getReadsFromCpGs(int min_cpg_coverage, List<RefCpG> list) {
		Set<MappedRead> mappedReadSet = new HashSet<>();
		list.forEach((refCpG) -> {
			refCpG.getCpGList().forEach((CpG cpg) -> mappedReadSet.add(cpg.getMappedRead()));
			assert refCpG.getCpGList().size() >= min_cpg_coverage; // make sure coverage is correct
		});
		return mappedReadSet;
	}

	/**
	 * write single interval output in mapped read format like the original data *
	 */
	private void writeMappedReadInInterval(String intervalFolderName, RefChr refChr, int startPos, int endPos,
	                                       Collection<MappedRead> mappedReadSet) throws IOException {
		BufferedWriter mappedReadWriter = new BufferedWriter(
				new FileWriter(String.format("%s/%s-%d-%d%s", intervalFolderName, refChr.getChr(), startPos, endPos,
						Constant.MAPPEDREADS_EXTENSION)));
		mappedReadWriter.write(String.format("ref:\t%s\n",
				refChr.getRefString().substring(startPos - INIT_POS, endPos + 1 - INIT_POS)));
		for (MappedRead mappedRead : mappedReadSet) {
			mappedReadWriter.write(mappedRead.outputString(startPos, endPos));
		}
		mappedReadWriter.close();
	}
}
