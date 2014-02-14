package ASM.DataType;

public class MappedRead {
	private String chrom;
	private String strand;
	private long start;
	private long end;
	private String sequence;
	private long id;

	public MappedRead(String chrom, String strand, long start, long end, String sequence, long id) {
		super();
		this.chrom = chrom;
		this.strand = strand;
		this.start = start;
		this.end = end;
		this.sequence = sequence;
		this.id = id;
	}

	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public String getStrand() {
		return strand;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}

	public long getStart() {
		return start;
	}

	public void setStart(long start) {
		this.start = start;
	}

	public long getEnd() {
		return end;
	}

	public void setEnd(long end) {
		this.end = end;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

	@Override
	public String toString() {
		return "MappedRead [chrom=" + chrom + ", strand=" + strand + ", start=" + start + ", end=" + end + ", sequence=" + sequence + ", id=" + id
				+ "]";
	}

	public String toRead() {
		return String.format("%s\t%s\t%d\t%d\t%s\t%d\n", chrom, strand, start, end, sequence, id);
	}
}
