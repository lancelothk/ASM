package edu.cwru.cbc.ASM.commons.Methylation;

/**
 * Created by lancelothk on 5/27/14.
 * Reference chromosome. Contains chr name and reference string.
 */
public class RefChr {
	private String chr;
	private String refString;

	/**
	 * @param refString should be upperCase string and without space in it.
	 */
	public RefChr(String chr, String refString) {
		this.chr = chr;
		this.refString = refString;
	}

	public String getChr() {
		return chr;
	}

	// The refString is in upper case.
	public String getRefString() {
		return refString;
	}

	public int getStart() {
		return 0;
	}

	public int getEnd() {
		return refString.length() - 1;
	}
}
