package edu.cwru.cbc.ASM.commons.DataType;

/**
 * Created by kehu on 11/12/14.
 */
public class EpiRead {
    private String chr;
    private int cpgOrder;
    private int cpgPos;
    private String cpgSeq;
    private String id;

    /**
     * Epiread format used in Amrfinder
     *
     * @param chr
     * @param cpgOrder
     * @param cpgSeq
     */
    public EpiRead(String chr, int cpgOrder, String cpgSeq) {
        this.chr = chr;
        this.cpgOrder = cpgOrder;
        this.cpgSeq = cpgSeq;
    }

    /**
     * Extended Epiread format
     *
     * @param chr
     * @param cpgOrder
     * @param cpgPos
     * @param cpgSeq
     * @param id
     */
    public EpiRead(String chr, int cpgOrder, int cpgPos, String cpgSeq, String id) {
        this.chr = chr;
        this.cpgOrder = cpgOrder;
        this.cpgPos = cpgPos;
        this.cpgSeq = cpgSeq;
        this.id = id;
    }

    public String epireadFormat() {
        if (!cpgSeq.contains("C") && !cpgSeq.contains("T")) {
            // if only contains N
            return "";
        } else {
            // use StringBuilder here is much faster than String.format
            StringBuilder sb = new StringBuilder();
            sb.append(chr);
            sb.append("\t");
            sb.append(cpgOrder);
            sb.append("\t");
            sb.append(cpgSeq);
            sb.append("\n");
            return sb.toString();
        }
    }

    public String extEpireadFormat() {
        if (!cpgSeq.contains("C") && !cpgSeq.contains("T")) {
            // if only contains N
            return "";
        } else {
            // use StringBuilder here is much faster than String.format
            StringBuilder sb = new StringBuilder();
            sb.append(chr);
            sb.append("\t");
            sb.append(cpgOrder);
            sb.append("\t");
            sb.append(cpgPos);
            sb.append("\t");
            sb.append(cpgSeq);
            sb.append("\t");
            sb.append(id);
            sb.append("\n");
            return sb.toString();
        }
    }
}
