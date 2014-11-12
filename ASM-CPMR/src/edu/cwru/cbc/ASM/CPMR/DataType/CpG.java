package edu.cwru.cbc.ASM.CPMR.DataType;

/**
 * Created by lancelothk on 5/27/14.
 * Store methylation information of CpG site in a read
 */
public class CpG {
	private MappedRead mappedRead; // link to MappedRead
	private CpGSite cpGSite; // link to CpGSite
	private MethylStatus methylStatus;

	public CpG(MappedRead mappedRead, CpGSite cpGSite) {
		this.mappedRead = mappedRead;
		this.cpGSite = cpGSite;
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

	public CpGSite getCpGSite() {
		return cpGSite;
	}

	public int getPos(){
		return this.cpGSite.getPos();
	}
}
