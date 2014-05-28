package edu.cwru.cbc.ASM.split.DataType;

/**
 * Created by lancelothk on 5/27/14.
 * Store methylation information of CpG site in a read
 */
public class CpG {
	private MappedRead mappedRead; // link to MappedRead
	private CpGSite cpGSite; // link to CpGSite
	private boolean isMethylated;

	public CpG(MappedRead mappedRead, CpGSite cpGSite) {
		this.mappedRead = mappedRead;
		this.cpGSite = cpGSite;
	}

	public boolean isMethylated() {
		return isMethylated;
	}

	public void setMethylated(boolean isMethylated) {
		this.isMethylated = isMethylated;
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
