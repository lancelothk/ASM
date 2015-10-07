package edu.cwru.cbc.ASM.commons.methylation;

import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

/**
 * Created by lancelothk on 5/27/14.
 * Store methylation information of CpG site in a read
 */
public class CpG {
	private MappedRead mappedRead; // link to MappedRead
	private RefCpG refCpG; // link to CpGSite
	private MethylStatus methylStatus;
	private int pos;

	public CpG(MappedRead mappedRead, RefCpG refCpG, MethylStatus methylStatus) {
		this.methylStatus = methylStatus;
		this.mappedRead = mappedRead;
		this.refCpG = refCpG;
		this.pos = refCpG.getPos();
	}

	public MethylStatus getMethylStatus() {
		return methylStatus;
	}

	public void setMethylStatus(MethylStatus methylStatus) {
		this.methylStatus = methylStatus;
	}

	public MappedRead getMappedRead() {
		return mappedRead;
	}

	public RefCpG getRefCpG() {
		return refCpG;
	}

	public int getPos() {
		return this.pos;
	}
}
