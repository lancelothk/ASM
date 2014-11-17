package edu.cwru.cbc.ASM.commons.DataType;

import com.google.common.base.Strings;

import java.util.ArrayList;
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
    private EpiRead epiRead;
    private int methylPolarity;

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
            case '-': {
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

    public String toString(int initialPos) throws IllegalArgumentException {
        if (this.start - initialPos > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("offset exceed max int!");
        }
        int offset = this.start - initialPos;
        return String.format("%s\t%s\t%d\t%d\t%s\t%s", this.chr, this.strand, this.start, this.end,
                             Strings.padStart(this.strand == '+' ? this.sequence : getComplementarySequence(),
                                              offset + this.sequence.length(), '.'), this.id);
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
        switch (cpg.getMethylStatus()) {
            case C:
                methylPolarity++;
                break;
            case T:
                methylPolarity--;
                break;
            case N:
                break;
            default:
                throw new RuntimeException("Unknown MethylStatus type!");
        }
    }

    public String outputString() {
        return String.format("%s\t%s\t%d\t%d\t%s\t%s\n", chr, strand, start, end, sequence, id);
    }

    public EpiRead getEpiRead() {
        if (epiRead == null && firstCpG != null) {
            return new EpiRead(chr, firstCpG.getRefCpG().getOrder(), firstCpG.getPos(), getCpGSeq(), id);
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

    public int getCpGCount() {
        return cpgList.size();
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

    public int getMethylPolarity() {
        return methylPolarity;
    }
}
