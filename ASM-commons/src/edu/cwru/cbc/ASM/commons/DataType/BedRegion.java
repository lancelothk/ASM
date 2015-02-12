package edu.cwru.cbc.ASM.commons.DataType;

/**
 * Created by kehu on 2/12/15.
 * Region in bed format.
 */
public class BedRegion {
    private String chr;
    private long start;
    private long end;
    private String name;

    public BedRegion(String chr, long start, long end, String name) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
    }
}
