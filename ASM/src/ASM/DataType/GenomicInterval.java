package ASM.DataType;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

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
    private BufferedWriter intervalWriter;

	public GenomicInterval(int start, int end, int length, int maxCount, double avgCount) {
		this.start = start;
		this.end = end;
		this.maxCount = maxCount;
		this.avgCount = avgCount;
		this.length = length;
	}

    public void initializeWriter(String intervalFolderName, String chr) throws IOException {
        intervalWriter = new BufferedWriter(
                new FileWriter(String.format("%s/%s-%d-%d.reads", intervalFolderName, chr, start, end)));
    }

    public void write(String toWrite) throws IOException {
        if (intervalWriter != null) {
            intervalWriter.write(toWrite);
        }
    }

    public void closeWriter() throws IOException {
        intervalWriter.close();
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
