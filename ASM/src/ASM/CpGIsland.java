package ASM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class CpGIsland {
	private String chrom;
	private long start;
	private long end;
	private String name;
	private int length;
	private int cpgNum;
	private int gcNum;
	private double perCpg;
	private double perGc;
	private double obsExp;
	private BufferedWriter bufferedWriter;

	public CpGIsland(String chrom, long start, long end, String name, int length, int cpgNum, int gcNum, double perCpg, double perGc, double obsExp) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.name = name;
		this.length = length;
		this.cpgNum = cpgNum;
		this.gcNum = gcNum;
		this.perCpg = perCpg;
		this.perGc = perGc;
		this.obsExp = obsExp;
	}

	public void initWriter(String outputPath) throws IOException {
		this.bufferedWriter = new BufferedWriter(new FileWriter(outputPath + "/" + this.name, true));
	}

	public void closeWriter() throws IOException {
		if (bufferedWriter != null) {
			this.bufferedWriter.close();
		}
	}

	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public long getStart() {
		return start;
	}

	public void setStart(long start) {
		this.start = start;
	}

	public long getEnd() {
		return end;
	}

	public void setEnd(long end) {
		this.end = end;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public int getCpgNum() {
		return cpgNum;
	}

	public void setCpgNum(int cpgNum) {
		this.cpgNum = cpgNum;
	}

	public int getGcNum() {
		return gcNum;
	}

	public void setGcNum(int gcNum) {
		this.gcNum = gcNum;
	}

	public double getPerCpg() {
		return perCpg;
	}

	public void setPerCpg(double perCpg) {
		this.perCpg = perCpg;
	}

	public double getPerGc() {
		return perGc;
	}

	public void setPerGc(double perGc) {
		this.perGc = perGc;
	}

	public double getObsExp() {
		return obsExp;
	}

	public void setObsExp(double obsExp) {
		this.obsExp = obsExp;
	}

	public void write(String str) throws IOException{
		bufferedWriter.write(str);
	}

	@Override
	public String toString() {
		return "CpGIsland [chrom=" + chrom + ", start=" + start + ", end=" + end + ", name=" + name + ", length=" + length + ", cpgNum=" + cpgNum
				+ ", gcNum=" + gcNum + ", perCpg=" + perCpg + ", perGc=" + perGc + ", obsExp=" + obsExp + "]";
	}

}
