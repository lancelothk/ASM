package ASM.DataType;

import com.google.common.base.Strings;

public class Read {
	private String chr;
	private char strand;
	private long start;
	private long end;
	private String content;
	private String id;

	public Read(String chr, char strand, long start, long end, String content, String id) {
		super();
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.end = end;
        this.content = content.toUpperCase();
        this.id = id;
	}

	public String getChr() {
		return chr;
	}

	public char getStrand() {
		return strand;
	}

	public long getStart() {
		return start;
	}

	public long getEnd() {
		return end;
	}

	public String getContent() {
		return content;
	}

	public String getId() {
		return id;
	}

	public String getComplementaryContent() {
		StringBuilder stringBuilder = new StringBuilder();
		for (char c : content.toCharArray()) {
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
				default:
					break;
			}
			stringBuilder.append(c);
		}
		return stringBuilder.toString();
	}

	public String toString(long initialPos) throws IllegalArgumentException {
		// TODO replace ledtPad with Guava pad
		if (this.start - initialPos > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("offset exceed max int!");
		}
		int offset = (int) (this.start - initialPos);
//        return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.end,
//                             Strings.padStart(this.content, offset + this.content.length(), '.'), this.id);
		return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.end, Strings.padStart(this.strand == '+' ? this.content : getComplementaryContent(), offset + this.content.length(), '.'), this.id);
    }
}
