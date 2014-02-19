package ASM.Execution;

import ASM.DataType.MappedRead;

import java.io.*;
import java.util.BitSet;

/**
 * Created by lancelothk on 2/17/14.
 */
public class CheckInterval {
	public static final int CHR6SIZE = 170899992;

	public static void main(String[] args) throws IOException {
		BitSet chrBitSet = new BitSet(CHR6SIZE);
		String chr = "chr6";
		String mappedReadFileName = "/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/reads_bs_i90_r1.mapped_chr6";
		String outputFileName = "/home/ke/test/checkInterval/" + chr;
		checkInterval(chrBitSet, chr, mappedReadFileName, outputFileName);
	}

	public static void checkInterval(BitSet chrBitSet, String chr, String mappedReadFileName, String outputFileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(mappedReadFileName));
		String line;
		String[] items;
		long counter = 0;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split("\t");
			if (items[0].startsWith("chrom")) {
				continue;
			}
			MappedRead mappedRead = new MappedRead(items[0], items[1], Long.parseLong(items[2]), Long.parseLong(items[3]), items[4],
					Long.parseLong(items[5]));
			chrBitSet.set((int)mappedRead.getStart(), (int)mappedRead.getEnd());
			counter++;
			if (counter % 1000000 == 0) {
				System.out.println(counter);
			}
		}
		bufferedReader.close();

		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		int clearIndex = 0, setIndex = 0;
		while (clearIndex <= chrBitSet.size() && setIndex <= chrBitSet.size()) {
			setIndex = chrBitSet.nextSetBit(clearIndex);
			if (setIndex < 0){
				break;
			}
			clearIndex = chrBitSet.nextClearBit(setIndex);
			bufferedWriter.write(String.format("%s\t%d\t%d\t%d\n", chr, setIndex, clearIndex - 1, clearIndex - setIndex));
		}
		bufferedWriter.close();
	}
}
