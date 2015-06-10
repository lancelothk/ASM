package edu.cwru.cbc.ASM.commons.GenomicInterval;

import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import org.apache.commons.lang3.math.NumberUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/12/15.
 * Genomic Region in Bed format.
 */
public class BedInterval extends GenomicIntervalBase implements Comparable<BedInterval> {
	private String name;
	private List<RefCpG> refCpGList;
	private boolean[] allelePattern;
	private boolean isPositive;
	private Set<BedInterval> intersectedRegions;

	public BedInterval(String chr, int start, int end, String name, boolean isPositive) {
		this(chr, start, end, name);
		this.isPositive = isPositive;
	}

	public BedInterval(String chr, int start, int end, String name) {
		super(chr, start, end);
		if (!chr.startsWith("chr") && NumberUtils.isNumber(chr)) {
			this.chr = "chr" + chr;
		}
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

	public void addIntersectedRegion(BedInterval region) {
		intersectedRegions.add(region);
	}

	public Set<BedInterval> getIntersectedRegions() {
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


	public String toBedWithIntersectionString() {
		return String.format("%s\t%d\t%d\t%s\t%s\t%b", this.chr, this.start, this.end, this.name, uniqueRegionnames(intersectedRegions), isIntersected());
	}

	private String uniqueRegionnames(Set<BedInterval> intersectedRegions) {
		Set<String> intersectionRegionNameSet = intersectedRegions.stream().map(BedInterval::getName).collect(
				Collectors.toSet());
		StringBuilder sb = new StringBuilder();
		for (String s : intersectionRegionNameSet) {
			sb.append(s).append(',');
		}
		return sb.toString();
	}

	@Override
	public int compareTo(BedInterval o) {
		int startDiff = this.getStart() - o.getStart();
		int endDiff = this.getEnd() - o.getEnd();
		return startDiff == 0 ? endDiff : startDiff;
	}
}
