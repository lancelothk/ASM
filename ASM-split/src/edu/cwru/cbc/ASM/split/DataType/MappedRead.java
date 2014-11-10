package edu.cwru.cbc.ASM.split.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/27/14.
 * Store mapped read in binary form to represent methylation status.
 */
public class MappedRead {
	private String chr;
	private String strand;
	private String id;
	private String sequence;
	private int start;
	private int end;
	private List<CpG> cpgList;

	public MappedRead(String chr, String strand, String id) {
		this.chr = chr;
		this.strand = strand;
		this.id = id;
		this.cpgList = new ArrayList<>();
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public String getChr() {
		return chr;
	}

	public String getStrand() {
		return strand;
	}

	public String getId() {
		return id;
	}

//	public int getStart() {
//		return cpgList.stream().min((a, b)-> a.getPos() - b.getPos()).get().getPos();
//	}
//
//	public int getEnd() {
//		return cpgList.stream().max((a, b)-> a.getPos() - b.getPos()).get().getPos();
//	}

	public void addCpG(CpG cpg){
		this.cpgList.add(cpg);
	}

	public int getCpGCount(){
		return cpgList.size();
	}

	public List<CpG> getCpgList() {
		return cpgList;
	}

	public String toWriteString(){
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start, end, sequence, id);
	}
}
