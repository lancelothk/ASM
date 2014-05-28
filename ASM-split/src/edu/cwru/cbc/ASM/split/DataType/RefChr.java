package edu.cwru.cbc.ASM.split.DataType;

/**
 * Created by lancelothk on 5/27/14.
 * Reference chromosome. Contains chr name and reference string.
 */
public class RefChr {
	private String chr;
	private String ref;

	public RefChr(String chr, String refString) {
		this.chr = chr;
		this.ref = refString;
	}

	public String getChr() {
		return chr;
	}

	public String getRefString() {
		return ref;
	}
}
