package edu.cwru.cbc.ASM.CPMR.DataType;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;

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
    private EpiRead epiRead;

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

    public EpiRead getEpiRead() {
        if (epiRead == null && firstCpG != null) {
            return new EpiRead(chr, firstCpG.getCpGSite().getId(), firstCpG.getPos(), getCpGSeq(), id);
        } else {
            return epiRead;
        }
    }

    private String getCpGSeq() {
        StringBuilder sb = new StringBuilder();
        for (CpG cpG : cpgList) {
            sb.append(cpG.getMethylStatus());
        }
        return sb.toString();
    }
}
