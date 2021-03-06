package edu.cwru.cbc.ASM.commons.methylation;

import javax.annotation.Nonnull;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

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
	private double p_value = -1;

	public RefCpG(int pos) {
		this.pos = pos;
		this.cpGList = new ArrayList<>();
	}

	public double getP_value() {
		return p_value;
	}

	public void setP_value(double p_value) {
		this.p_value = p_value;
	}

	/**
	 * check if parameter RefCpG have shared reads  with current RefCpG *
	 */
	public boolean hasCommonRead(RefCpG next, int min_cpg_coverage) {
		int count = 0;
		for (CpG cCpg : this.getCpGList()) {
			for (CpG nCpg : next.getCpGList()) {
				if (cCpg.getMappedRead() == nCpg.getMappedRead()) {
					count++;
					if (count >= min_cpg_coverage) {
						return true;
					}
				}
			}
		}
		return false;
	}

	public List<CpG> getCpGList() {
		return cpGList;
	}

	public void assignIndex(int order) {
		this.index = order;
	}

	public int getIndex() {
		return index;
	}

	public void addCpG(CpG cpg) {
		this.cpGList.add(cpg);
		if (cpg.getMethylStatus() == MethylStatus.C) {
			methylCount++;
			coveredCount++;
		} else {
			coveredCount++;
		}
	}

	public int getCpGCoverage() {
		return cpGList.size();
	}

	public boolean hasPartialMethyl(double partial_methyl_threshold) {
		double m = 0, n = 0;
		for (CpG cpG : cpGList) {
			if (cpG.getMethylStatus() == MethylStatus.C) {
				m++;
			} else if (cpG.getMethylStatus() == MethylStatus.T) {
				n++;
			}
		}
		double rate = m > n ? n / (m + n) : m / (m + n);
		return rate >= partial_methyl_threshold;
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
	public int compareTo(@Nonnull RefCpG o) {
		return this.getPos() - o.getPos();
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof RefCpG)) return false;
		RefCpG refCpG = (RefCpG) o;
		return Objects.equals(pos, refCpG.pos);
	}

	@Override
	public int hashCode() {
		return Objects.hash(pos);
	}

	public int getPos() {
		return pos;
	}
}
