package edu.cwru.cbc.ASM.detect.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 3/5/14.
 * Genome Interval Data Structure
 * Intervals should not overlap with each other.
 */
public class GenomicInterval {
	private int start;
	private int end;
	private int maxCount;
	private double avgCount;
	private int readCount;
	private int length;
    private List<String> readList;

	public GenomicInterval(int start, int end, int length, int maxCount, double avgCount) {
		this.start = start;
		this.end = end;
		this.maxCount = maxCount;
		this.avgCount = avgCount;
		this.length = length;
        this.readList = new ArrayList<>();
    }

    public void addRead(String read) {
        this.readList.add(read);
    }

    public List<String> getReadList() {
        return readList;
    }

    public void setReadCount(int readCount) {
        this.readCount = readCount;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public int getReadCount() {
		return readCount;
	}

	public void incrementReadCount() {
		this.readCount++;
	}

	public int getMaxCount() {
		return maxCount;
	}

	public void setMaxCount(int maxCount) {
		this.maxCount = maxCount;
	}

	public double getAvgCount() {
		return avgCount;
	}

	public void setAvgCount(double avgCount) {
		this.avgCount = avgCount;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}
}
