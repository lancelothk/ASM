package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * Store reference CpG site information in reference genome.
 */
public class RefCpG implements Comparable<RefCpG> {
	private int pos;
	private int index;
	private List<CpG> cpGList;
	private int methylCount;
	private int coveredCount;

	public RefCpG(int pos) {
		this.pos = pos;
		this.cpGList = new ArrayList<>();
	}

	public RefCpG(int pos, MethylStatus methylStatus) {
		this.pos = pos;
		this.cpGList = new ArrayList<>();
		if (methylStatus == MethylStatus.C) {
			methylCount++;
			coveredCount++;
		} else {
			coveredCount++;
		}
	}

	/**
	 * check if parameter RefCpG have shared reads  with current RefCpG *
	 */
	public boolean hasCommonRead(RefCpG next) {
		for (CpG cCpg : this.getCpGList()) {
			for (CpG nCpg : next.getCpGList()) {
				if (cCpg.getMappedRead() == nCpg.getMappedRead()) {
					return true;
				}
			}
		}
		return false;
	}

	public int getPos() {
		return pos;
	}

	public void assignIndex(int order) {
		this.index = order;
	}

	public int getIndex() {
		return index;
	}

	public void addCpG(CpG cpg) {
		this.cpGList.add(cpg);
	}

	public int getCpGCoverage() {
		return cpGList.size();
	}

	public boolean hasPartialMethyl() {
		int m = 0, n = 0;
		for (CpG cpG : cpGList) {
			if (cpG.getMethylStatus() == MethylStatus.C) {
				m++;
			} else if (cpG.getMethylStatus() == MethylStatus.T) {
				n++;
			}
		}
		return m != 0 && n != 0;
	}

	public List<CpG> getCpGList() {
		return cpGList;
	}

	public void addMethylCount(int count) {
		methylCount += count;
		coveredCount += count;
	}

	public void addNonMethylCount(int count) {
		coveredCount += count;
	}

	public int getCoveredCount() {
		return coveredCount;
	}

	public double getMethylLevel() {
		return (double) methylCount / coveredCount;
	}

	public int getMethylCount() {
		return methylCount;
	}

	public int getNonMethylCount() {
		return coveredCount - methylCount;
	}

	public MethylStatus getMajorMethylStatus() {
		if (methylCount > coveredCount / 2.0) {
			return MethylStatus.C;
		} else if (methylCount < coveredCount / 2.0) {
			return MethylStatus.T;
		} else {
			// if methylCount == nonMethylCount, return unknown
			return MethylStatus.E;
		}
	}

	@Override
	public int compareTo(RefCpG o) {
		return this.getPos() - o.getPos();
	}
}
