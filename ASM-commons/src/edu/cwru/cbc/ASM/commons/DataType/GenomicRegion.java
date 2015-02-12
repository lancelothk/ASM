package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by kehu on 2/12/15.
 * Genomic Region in bed format.
 */
public class GenomicRegion {
    private String chr;
    private long start;
    private long end;
    private String name;
    private List<RefCpG> refCpGList;
    private boolean[] allelePattern;


    public GenomicRegion(String chr, long start, long end, String name) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.refCpGList = new ArrayList<>();
    }

    public String getChr() {
        return chr;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    public String getName() {
        return name;
    }

    public void addRefCpG(RefCpG refCpG) {
        this.refCpGList.add(refCpG);
    }

    public List<RefCpG> getRefCpGList() {
        return refCpGList;
    }

    public void generateRandomAllelePattern() {
        this.allelePattern = new boolean[refCpGList.size()];
        Random rand = new Random();
        for (int i = 0; i < allelePattern.length; i++) {
            allelePattern[i] = rand.nextBoolean();
        }

    }
}
