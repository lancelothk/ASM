package edu.cwru.cbc.ASM.CPMR.DataType;

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
    private CpG firstCpG;

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

    public void addCpG(CpG cpg) {
        this.cpgList.add(cpg);
        if (firstCpG == null || cpg.getPos() < firstCpG.getPos()) {
            firstCpG = cpg;
        }
    }

    public int getCpGCount() {
        return cpgList.size();
    }

    public List<CpG> getCpgList() {
        return cpgList;
    }

    public String outputString() {
        return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start, end, sequence, id);
    }

    private String getCpGSeq() {
        StringBuilder sb = new StringBuilder();
        for (CpG cpG : cpgList) {
            sb.append(cpG.getMethylStatus());
        }
        return sb.toString();
    }

    public String extEpireadFormat() {
        if (firstCpG == null) {
            return "";
        } else {
            // use StringBuilder here is much faster than String.format
            StringBuilder sb = new StringBuilder();
            sb.append(chr);
            sb.append("\t");
            sb.append(firstCpG.getCpGSite().getId());
            sb.append("\t");
            sb.append(firstCpG.getPos());
            sb.append("\t");
            sb.append(getCpGSeq());
            sb.append("\t");
            sb.append(id);
            sb.append("\n");
            return sb.toString();
//            return String.format("%s\t%d\t%d\t%s\t%s\n", chr, firstCpG.getCpGSite().getId(), firstCpG.getPos(), getCpGSeq(), id);
        }
    }
}
