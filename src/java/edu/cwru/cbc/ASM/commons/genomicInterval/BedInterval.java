package edu.cwru.cbc.ASM.commons.genomicInterval;

import edu.cwru.cbc.ASM.commons.methylation.RefCpG;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/12/15.
 * Genomic Region in Bed format.
 */
public class BedInterval extends GenomicIntervalBase {
	private String name;
	private List<RefCpG> refCpGList;
	private boolean[] allelePattern;
	private Set<BedInterval> intersectedRegions;

	public BedInterval(String chr, int start, int end, String name) {
		super(chr, start, end);
		this.name = name;
		this.refCpGList = new ArrayList<>();
		this.intersectedRegions = new HashSet<>();
	}

	public void addIntersectedRegion(BedInterval region) {
		intersectedRegions.add(region);
	}

	public Set<BedInterval> getIntersectedRegions() {
		return intersectedRegions;
	}

	public boolean isIntersected() {
		return intersectedRegions.size() > 0;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
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


	public String toBedWithIntersectionString() {
		return String.format("%s\t%d\t%d\t%s\t%s\t%b", this.chr, this.start, this.end, this.name, uniqueRegionNames(intersectedRegions), isIntersected());
	}

	private String uniqueRegionNames(Set<BedInterval> intersectedRegions) {
		Set<String> intersectionRegionNameSet = intersectedRegions.stream().map(BedInterval::getName).collect(
				Collectors.toSet());
		StringBuilder sb = new StringBuilder();
		for (String s : intersectionRegionNameSet) {
			sb.append(s).append(',');
		}
		return sb.toString();
	}
}
