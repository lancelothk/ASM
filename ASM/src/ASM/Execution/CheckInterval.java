package ASM.Execution;

import ASM.DataType.MappedRead;
import ASM.Utils.IntervalChekingLineProcessor;
import ASM.Utils.MappedReadFileLineProcessor;
import com.google.common.io.CharStreams;

import java.io.*;
import java.util.BitSet;

/**
 * Created by lancelothk on 2/17/14.
 */
public class CheckInterval {
	public static final int CHR6SIZE = 170899992;

	public static void main(String[] args) throws IOException {
		String chr = "chr6";
		String mappedReadFileName = "/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/reads_bs_i90_r1.mapped_chr6";
		String outputFileName = "/home/ke/test/checkInterval/" + chr;
		checkInterval(CHR6SIZE, chr, mappedReadFileName, outputFileName);
	}

	public static void checkInterval(int chrBitSize, String chr, String mappedReadFileName, String outputFileName) throws IOException {
		BitSet chrBitSet = CharStreams.readLines(new BufferedReader((new FileReader(mappedReadFileName))), new IntervalChekingLineProcessor(chrBitSize));
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
