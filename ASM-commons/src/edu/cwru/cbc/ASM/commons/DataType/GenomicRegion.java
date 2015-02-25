package edu.cwru.cbc.ASM.commons.DataType;

import java.util.ArrayList;
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


    public GenomicRegion(String chr, int start, int end, String name) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.refCpGList = new ArrayList<>();
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

	public String toBedString(){
		return String.format("%s\t%d\t%d\t%s", this.chr, this.start, this.end, this.name);
	}

    @Override
    public int compareTo(GenomicRegion o) {
        int startDiff = this.getStart() - o.getStart();
        int endDiff = this.getEnd() - o.getEnd();
        return startDiff == 0 ? endDiff : startDiff;
    }
}
