
public class Read {
	private String chr;
	private String strand;
	private long start;
	private long end;
	private String content;
	public Read(String chr, String strand, long start, long end, String content) {
		super();
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.end = end;
		this.content = content;
	}
	public String getChr() {
		return chr;
	}
	public String getStrand() {
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
}
