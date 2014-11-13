package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * Store reference CpG site information in reference genome.
 */
public class RefCpG {
	private int pos;
	private int id;
	private List<CpG> cpGList;

	public RefCpG(int pos) {
		this.pos = pos;
		this.cpGList = new ArrayList<>();
	}

	public int getPos() {
		return pos;
	}

	public void assignId(int id) {
		this.id = id;
	}

	public int getId() {
		return id;
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
}
