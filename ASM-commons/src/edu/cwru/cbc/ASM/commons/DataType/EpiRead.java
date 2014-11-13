package edu.cwru.cbc.ASM.commons.DataType;

/**
 * Created by kehu on 11/12/14.
 * Should be used for single end data. If pair-end, cpgSeq may contain gaps.
 */
public class EpiRead {
    private String chr;
    private int cpgOrder;
    private int cpgPos;
    private String cpgSeq;
    private String id;

    /**
     * Epiread format used in Amrfinder
     */
    public EpiRead(String chr, int cpgOrder, String cpgSeq) {
        this.chr = chr;
        this.cpgOrder = cpgOrder;
        this.cpgSeq = cpgSeq;
    }

    /**
     * Extended Epiread format
     */
    public EpiRead(String chr, int cpgOrder, int cpgPos, String cpgSeq, String id) {
        this.chr = chr;
        this.cpgOrder = cpgOrder;
        this.cpgPos = cpgPos;
        this.cpgSeq = cpgSeq;
        this.id = id;
    }

    /**
     * Read epiread format from raw line
     * @param line
     * @param format
     */
    public EpiRead(String line, EpiReadFormat format) {
        String[] items = line.split("\t");
        switch (format) {
            case epiread:
                if (items.length != 3) {
                    throw new RuntimeException("invalid line for epiread format!");
                }
                new EpiRead(items[0], Integer.parseInt(items[1]), items[3]);
                break;
            case extEpiread:
                if (items.length != 5) {
                    throw new RuntimeException("invalid line for extEpiread format!");
                }
                new EpiRead(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]), items[3], items[4]);
                break;
            default:
                throw new RuntimeException("Unknown epiread format!");
        }
    }

    public String getChr() {
        return chr;
    }

    public int getCpgOrder() {
        return cpgOrder;
    }

    public int getCpgPos() {
        return cpgPos;
    }

    public String getCpgSeq() {
        return cpgSeq;
    }

    public String getId() {
        return id;
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

    public enum EpiReadFormat {
        epiread, extEpiread
    }
}
