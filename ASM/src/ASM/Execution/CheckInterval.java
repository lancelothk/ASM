package ASM.Execution;

import ASM.DataType.ChrCoverageSummary;
import ASM.Utils.IntervalChekingLineProcessor;
import com.google.common.io.CharStreams;

import java.io.*;

/**
 * Created by lancelothk on 2/17/14.
 * Check mapped reads interval in Chromosome
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
		ChrCoverageSummary chrCoverageSummary = CharStreams.readLines(new BufferedReader((new FileReader(mappedReadFileName))), new IntervalChekingLineProcessor(chrBitSize));
		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		chrCoverageSummary.writeIntervalCoverageSummary(chr, bufferedWriter);
		bufferedWriter.close();
	}
}
