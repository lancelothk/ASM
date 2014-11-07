package edu.cwru.cbc.ASM.split.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * Store CpG site information in reference genome.
 */
public class CpGSite {
	private int pos;
	private int cpgSiteID;
	private List<CpG> cpGList;

	public CpGSite(int pos, int cpgSiteID) {
		this.pos = pos;
		this.cpgSiteID = cpgSiteID;
		this.cpGList = new ArrayList<>();
	}

	public int getPos() {
		return pos;
	}

	public int getCpgSiteID() {
		return cpgSiteID;
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
            if (cpG.isMethylated()) {
                m++;
            } else {
                n++;
            }
        }
        return m != 0 && n != 0;
    }

	public List<CpG> getCpGList() {
		return cpGList;
	}
}
