package ASM.DataType;

import java.util.ArrayList;
import java.util.List;

public class MappedRead {
	private String chrom;
	private String strand;
	private long start;
	private long end;
	private String sequence;
	private long id;
	private List<CpGSite> cpgList;
    private CpGSite firstCpG;
    private CpGSite lastCpG;

	public MappedRead(String chrom, String strand, long start, long end, String sequence, long id) {
		super();
		this.chrom = chrom;
		this.strand = strand;
		this.start = start;
		this.end = end;
        this.sequence = sequence.toUpperCase();
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

	public List<CpGSite> getCpgList() {
		return cpgList;
	}

    public void addCpG(CpGSite cpg) {
        if (cpgList == null) {
            cpgList = new ArrayList<>();
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

    public CpGSite getFirstCpG() {
        return firstCpG;
    }

    public CpGSite getLastCpG() {
        return lastCpG;
    }

    public void setCpgList(List<CpGSite> cpgList) {
        this.cpgList = cpgList;
	}

	@Override
	public String toString() {
		return "MappedRead [chrom=" + chrom + ", strand=" + strand + ", start=" + start + ", end=" + end + ", sequence=" + sequence + ", id=" + id
				+ "]";
	}

	public String toRead() {
		return String.format("%s\t%s\t%d\t%d\t%s\t%d\n", chrom, strand, start, end, sequence, id);
	}

    public boolean getCpGMethylStatus(long pos) {
        if (strand.equals("+")) {
            // CG is methylated, TG is non-methylated, others are ignored.
            char bp = sequence.charAt((int) (pos - start));
            if (bp == 'C') {
                return true;
            } else if (bp == 'T') {
                return false;
            } else {
                // TODO should ignore non C/T case in CpG list
                return false;
            }
        } else if (strand.equals("-")) {
            // GC is methylated, GT is non-methylated, others are ignored.
            char bp = sequence.charAt((int) (pos - start + 1));
            if (bp == 'G') {
                return true;
            } else if (bp == 'A') {
                return false;
            } else {
                // TODO should ignore non C/T case in CpG list
                return false;
            }
        } else {
            throw new RuntimeException("illegal strand!");
        }
    }
}
