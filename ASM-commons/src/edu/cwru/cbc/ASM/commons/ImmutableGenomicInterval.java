package edu.cwru.cbc.ASM.commons;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;

import java.util.List;

/**
 * Created by lancelothk on 6/8/15.
 * Store immutable genomic interval for passing and output.
 * start pos is 0-based, end pos is 0-based.
 */
public class ImmutableGenomicInterval extends ImmutableGenomicIntervalBase {

	private final String refString;
	private ImmutableList<RefCpG> refCpGList;
	private ImmutableList<MappedRead> mappedReadList;

	public ImmutableGenomicInterval(String chr, String refString, int start, int end, List<RefCpG> refCpGList,
	                                List<MappedRead> mappedReadList) {
		super(chr, start, end);
		this.refString = refString;
		this.refCpGList = ImmutableList.copyOf(refCpGList);
		this.mappedReadList = ImmutableList.copyOf(mappedReadList);
	}

	public String getRefString() {
		return refString;
	}

	public ImmutableList<RefCpG> getRefCpGList() {
		return refCpGList;
	}

	public ImmutableList<MappedRead> getMappedReadList() {
		return mappedReadList;
	}
}
