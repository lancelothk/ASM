package edu.cwru.cbc.ASM.commons.sequence;

import com.google.common.base.Strings;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.function.BiConsumer;

/**
 * Created by lancelothk on 11/12/14.
 * Used to represent mapped read.
 */
public class MappedRead implements Comparable<MappedRead> {
	private final String id;
	private final String sequence;
	private final String chr;
	private final char strand;
	private int start;
	private List<CpG> cpgList;

	public MappedRead(String chr, char strand, int start, String sequence, String id) {
		this.id = id;
		this.sequence = sequence;
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.cpgList = new ArrayList<>();
	}

	public static void CpGInRead(char strand, String sequence, int startPos, HashIntObjMap<RefCpG> refMap,
	                             BiConsumer<RefCpG, MethylStatus> action) {
		int start = strand == '+' ? startPos : startPos - 1;
		int end = strand == '+' ? startPos + sequence.length() - 1 : startPos + sequence.length() - 2;
		for (int i = start; i <= end; i++) {
			if (refMap.containsKey(i)) {
				MethylStatus methylStatus = getMethylStatus(i, strand, sequence, startPos);
				RefCpG refCpG = refMap.get(i);
				if (refCpG != null && methylStatus != MethylStatus.N) {
					action.accept(refCpG, methylStatus);
					i++;
				}
			}
		}
	}

	public static MethylStatus getMethylStatus(int pos, char strand, String sequence, int start) {
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

	public static String getComplementarySequence(String sequence) {
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
				case '-':
					c = '-';
					break;
				default:
					throw new RuntimeException("invalid character in sequence!");
			}
			stringBuilder.append(c);
		}
		return stringBuilder.toString();
	}

	/**
	 * This method associate RefCpG <-> CpG <-> MappedRead
	 */
	public int generateCpGsInRead(HashIntObjMap<RefCpG> refMap) {
		CpGInRead(strand, sequence, start, refMap, (refCpG, methylStatus) -> {
			CpG cpg = new CpG(this, refCpG, methylStatus);
			this.cpgList.add(cpg);
			refCpG.addCpG(cpg);
		});
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

	public String toVisualizationString(int initialPos) throws IllegalArgumentException {
		int offset = this.start - initialPos;
		return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.getEnd(),
				Strings.padStart(this.strand == '+' ? this.sequence : getComplementarySequence(),
						offset + this.sequence.length(), '.'), this.id);
	}

	public String getComplementarySequence() {
		return getComplementarySequence(sequence);
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
		return this.cpgList.get(0);
	}

	public CpG getLastCpG() {
		return this.cpgList.get(this.cpgList.size() - 1);
	}

	public char getStrand() {
		return strand;
	}

	public String getId() {
		return id;
	}

	public String getSequence() {
		return sequence;
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

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof MappedRead)) return false;
		MappedRead that = (MappedRead) o;
		return Objects.equals(id, that.id);
	}

	@Override
	public int hashCode() {
		return Objects.hash(id);
	}
}
