package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by kehu on 2/12/15.
 * Genomic Region in bed format.
 */
public class GenomicRegion implements Comparable<GenomicRegion> {
    private String chr;
    private int start;
    private int end;
    private String name;
    private List<RefCpG> refCpGList;
    private boolean[] allelePattern;
    private boolean isIntersected;
    private boolean isPositive;

    public GenomicRegion(String chr, int start, int end, String name) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.refCpGList = new ArrayList<>();
    }

    public GenomicRegion(String chr, int start, int end, String name, boolean isPositive) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.isPositive = isPositive;
    }

    public boolean isPositive() {
        return isPositive;
    }

    public void setPositive(boolean isPositive) {
        this.isPositive = isPositive;
    }

    public boolean isIntersected() {
        return isIntersected;
    }

    public void setIntersected(boolean isIntersected) {
        this.isIntersected = isIntersected;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
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

    public boolean getRefMethylStatus(int index) {
        return allelePattern[index];
    }

    public String toPatternString() {
        return String.format("%d\t%d\t%s", start, end, Arrays.toString(allelePattern));
    }

    public String toBedString() {
        return String.format("%s\t%d\t%d\t%s", this.chr, this.start, this.end, this.name);
	}

    @Override
    public int compareTo(GenomicRegion o) {
        int startDiff = this.getStart() - o.getStart();
        int endDiff = this.getEnd() - o.getEnd();
        return startDiff == 0 ? endDiff : startDiff;
    }
}
