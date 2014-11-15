package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * Store reference CpG site information in reference genome.
 */
public class RefCpG implements Comparable<RefCpG> {
	private int pos;
	private int order;
	private List<CpG> cpGList;
	private int methylCount;
	private int coveredCount;

	public RefCpG(int pos) {
		this.pos = pos;
		this.cpGList = new ArrayList<>();
	}

	public int getPos() {
		return pos;
	}

	public void assignOrder(int order) {
		this.order = order;
	}

	public int getOrder() {
		return order;
	}

	public void addCpG(CpG cpg){
		this.cpGList.add(cpg);
	}

	public int getCoverage(){
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

	@Override
	public int compareTo(RefCpG o) {
		return this.getPos() - o.getPos();
	}
}
