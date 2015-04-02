package edu.cwru.cbc.ASM.commons;

import edu.cwru.cbc.ASM.commons.CpG.RefCpG;

import java.util.*;

/**
 * Created by kehu on 2/12/15.
 * Genomic Region in Bed format.
 */
public class GenomicInterval extends GenomicIntervalBase implements Comparable<GenomicInterval> {
	private String name;
	private List<RefCpG> refCpGList;
	private boolean[] allelePattern;
	private boolean isPositive;
	private Set<GenomicInterval> intersectedRegions;

	public GenomicInterval(String chr, int start, int end, String name, boolean isPositive) {
		this(chr, start, end, name);
		this.isPositive = isPositive;
	}

	public GenomicInterval(String chr, int start, int end, String name) {
		super(chr, start, end);
		this.name = name;
		this.refCpGList = new ArrayList<>();
		this.intersectedRegions = new HashSet<>();
	}

	public boolean isPositive() {
		return isPositive;
	}

	public void setPositive(boolean isPositive) {
		this.isPositive = isPositive;
	}

	public void addIntersectedRegion(GenomicInterval region) {
		intersectedRegions.add(region);
	}

	public Set<GenomicInterval> getIntersectedRegions() {
		return intersectedRegions;
	}

	public boolean isIntersected() {
		return intersectedRegions.size() != 0;
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

	public int length() {
		return this.getEnd() - this.getStart() + 1;
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
	public int compareTo(GenomicInterval o) {
		int startDiff = this.getStart() - o.getStart();
		int endDiff = this.getEnd() - o.getEnd();
		return startDiff == 0 ? endDiff : startDiff;
	}
}
