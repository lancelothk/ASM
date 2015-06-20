package edu.cwru.cbc.ASM.commons.sequence;

import com.google.common.base.Strings;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Created by lancelothk on 11/12/14.
 * Used to represent mapped read.
 */
public class MappedRead extends Sequence implements Comparable<MappedRead> {
	private final String chr;
	private final char strand;
	private int start;
	private List<CpG> cpgList;
	private CpG firstCpG;
	private CpG lastCpG;

	public MappedRead(String chr, char strand, int start, String sequence, String id) {
		super(id, sequence.toUpperCase());
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.cpgList = new ArrayList<>();
	}

	/**
	 * This method associate RefCpG <-> CpG <-> MappedRead
	 */
	public int generateCpGsInRead(Map<Integer, RefCpG> refMap) {
		int start = this.strand == '+' ? this.getStart() : this.getStart() - 1;
		int end = this.strand == '+' ? this.getEnd() : this.getEnd() - 1;
		for (int i = start; i <= end; i++) {
			if (refMap.containsKey(i)) {
				MethylStatus methylStatus = this.getMethylStatus(i);
				if (methylStatus != MethylStatus.N) {
					CpG cpg = new CpG(this, refMap.get(i), methylStatus);
					this.addCpG(cpg);
					refMap.get(i).addCpG(cpg);
					i++;
				}
			}
		}
		return this.getCpgList().size();
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	/**
	 * 0-based (inclusive the last base pair)
	 */
	public int getEnd() {
		return start + sequence.length() - 1;
	}

	private MethylStatus getMethylStatus(int pos) {
		// minus strand is the complementary string of plus strand.
		// only consider 'C' position in CpG. Ignore the char in 'G' position
		switch (strand) {
			case '+': {
				// CN is methylated, TN is non-methylated
				char cbp = sequence.charAt(pos - start);
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
		int offset = this.start - initialPos;
		return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.getEnd(),
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
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", this.chr, this.strand, this.start, this.getEnd(),
				new String(strArray), this.id);
	}

	// convert to MR format sequence
	public String toMRFormatString(int mismatch, char qualityChar) {
		char[] qualityStrArray = new char[this.sequence.length()];
		Arrays.fill(qualityStrArray, qualityChar);
		return String.format("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n", this.chr, this.start, this.getEnd(), this.id, mismatch,
				this.strand, this.sequence, new String(qualityStrArray));
	}

	private void addCpG(CpG cpg) {
		if (cpgList.size() == 0) {
			this.cpgList.add(cpg);
			firstCpG = cpg;
			lastCpG = cpg;
		} else {
			this.cpgList.add(cpg);
			if (cpg.getPos() > lastCpG.getPos()) {
				lastCpG = cpg;
			}
		}
	}

	public String toMappedReadString() {
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start, this.getEnd(), sequence, id);
	}

	public String toMappedReadString(int leftBound, int rightBound) {
		int leftOffset = start < leftBound ? leftBound - start : 0;
		int rightOffset = this.getEnd() > rightBound ? this.getEnd() - rightBound : 0;
		return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start + leftOffset, this.getEnd() - rightOffset,
				sequence.substring(leftOffset, sequence.length() - rightOffset), id);
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

	@Override
	public int compareTo(MappedRead other) {
		int startCompare = Integer.compare(this.getStart(), other.getStart());
		if (startCompare != 0) {
			return startCompare;
		}
		int endCompare = Integer.compare(this.getEnd(), other.getEnd());
		if (endCompare != 0) {
			return endCompare;
		}
		return this.getId().compareTo(other.getId());
	}
}
