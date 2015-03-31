package edu.cwru.cbc.ASM.commons.DataType;

import com.google.common.base.Strings;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by lancelothk on 11/12/14.
 * Used to represent mapped read.
 */
public class MappedRead {
	private String chr;
	private String sequence;
	private String id;
	private char strand;
	private int start;
	private int end;
	private List<CpG> cpgList;
	private CpG firstCpG;
	private CpG lastCpG;

	public MappedRead(String chr, char strand, int start, int end, String sequence, String id) {
		super();
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.end = end;
		this.sequence = sequence.toUpperCase();
		this.id = id;
		this.cpgList = new ArrayList<>();
	}


	public String getChr() {
		return chr;
	}

	public String getSequence() {
		return sequence;
	}

	public String getId() {
		return id;
	}

	public char getStrand() {
		return strand;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public MethylStatus getMethylStatus(int pos) {
		// only consider 'C' position in CpG. Ignore the char in 'G' position
		switch (strand) {
			case '+': {
				// CN is methylated, TN is non-methylated
				char cbp = sequence.charAt(pos - start);
				char gbp = sequence.charAt(pos - start + 1);
				if (cbp == 'C') {
					return MethylStatus.C;
				} else if (cbp == 'T') {
					return MethylStatus.T;
				} else {
					return MethylStatus.N;
				}
			}
			case '-': { // here the - strand is origianl one. Not the complementary one.
				// NC is methylated, NT is non-methylated
				char cbp = sequence.charAt(pos - start + 1);
				char gbp = sequence.charAt(pos - start);
				if (cbp == 'C') {  // for complementary bp, 'G'
					return MethylStatus.C;
				} else if (cbp == 'T') { // for complementary bp, 'A'
					return MethylStatus.T;
				} else {
					return MethylStatus.N;
				}
			}
			default:
				throw new RuntimeException("illegal strand!");
		}
	}

	public String toVisualizationString(int initialPos) throws IllegalArgumentException {
		if (this.start - initialPos > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("offset exceed max int!");
		}
		int offset = this.start - initialPos;
		return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.end,
				Strings.padStart(this.strand == '+' ? this.sequence : getComplementarySequence(),
						offset + this.sequence.length(), '.'), this.id);
	}

	public String getComplementarySequence() {
		StringBuilder stringBuilder = new StringBuilder();
		for (char c : sequence.toCharArray()) {
			switch (c) {
				case 'A':
					c = 'T';
					break;
				case 'T':
					c = 'A';
					break;
				case 'C':
					c = 'G';
					break;
				case 'G':
					c = 'C';
					break;
				case 'N':
					c = 'N';
					break;
				default:
					throw new RuntimeException("invalid character in sequence!");
			}
			stringBuilder.append(c);
		}
		return stringBuilder.toString();
	}

	// + strand     CG  TG
	// - strand     GC  AC  GT
	// - comp       CG  TG  CA
	public String toSimulationString() {
		char[] strArray = sequence.toCharArray();
		for (CpG cpG : cpgList) {
			if (cpG.getMethylStatus() == MethylStatus.C) {
				if (strand == '+') {
					strArray[cpG.getPos() - this.start] = 'C';
					strArray[cpG.getPos() - this.start + 1] = 'G';
				} else {
					strArray[cpG.getPos() - this.start + 1] = 'C';
					strArray[cpG.getPos() - this.start] = 'G';
				}
			}
			if (cpG.getMethylStatus() == MethylStatus.T) {
				if (strand == '+') {
					strArray[cpG.getPos() - this.start] = 'T';
					strArray[cpG.getPos() - this.start + 1] = 'G';
				} else {
					strArray[cpG.getPos() - this.start + 1] = 'T';
					strArray[cpG.getPos() - this.start] = 'G';
				}
			}
		}
		return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.end,
				new String(strArray), this.id);
	}

	// convert to MR format sequence
	public String toMRFormatString(int mismatch, char qualityChar) {
		char[] qualityStrArray = new char[this.sequence.length()];
		Arrays.fill(qualityStrArray, qualityChar);
		return String.format("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s", this.chr, this.start, this.end, this.id, mismatch,
				this.strand, this.sequence, new String(qualityStrArray));
	}

	public void addCpG(CpG cpg) {
		if (cpgList.size() == 0) {
			this.cpgList.add(cpg);
			firstCpG = cpg;
			lastCpG = cpg;
		} else {
			this.cpgList.add(cpg);
			if (cpg.getPos() < firstCpG.getPos()) {
				firstCpG = cpg;
			}
			if (cpg.getPos() > lastCpG.getPos()) {
				lastCpG = cpg;
			}
		}
	}

	public String outputString() {
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start, end, sequence, id);
	}

	public String outputString(int leftBound, int rightBound) {
		int leftOffset = start < leftBound ? leftBound - start : 0;
		int rightOffset = end > rightBound ? end - rightBound : 0;
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start + leftOffset, end - rightOffset,
				sequence.substring(leftOffset, sequence.length() - rightOffset), id);
	}

	private String getCpGSeq(int startCpGPos, int endCpGPos) {
		StringBuilder sb = new StringBuilder();
		for (CpG cpG : cpgList) {
			if (cpG.getPos() >= startCpGPos && cpG.getPos() <= endCpGPos) {
				sb.append(cpG.getMethylStatus());
			}
		}
		return sb.toString();
	}

	public List<CpG> getCpgList() {
		return cpgList;
	}

	public CpG getFirstCpG() {
		return firstCpG;
	}

	public CpG getLastCpG() {
		return lastCpG;
	}

}
