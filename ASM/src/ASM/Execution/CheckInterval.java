package ASM.Execution;

import ASM.DataType.CpGIsland;
import ASM.DataType.MappedRead;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;

/**
 * Created by lancelothk on 2/17/14.
 */
public class CheckInterval {
	public static final int CHR6SIZE = 170899992;

	public static void main(String[] args) {
		BitSet chrBitSet = new BitSet(CHR6SIZE);
		String mappedReadFileName = "";
	}

	public static void checkInterval(String mappedReadFileName, BitSet chrBitSet) throws IOException {

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

		int fromIndex = 1;
		int clearIndex = 0;
		while (clearIndex <= chrBitSet.size()) {
			//TODO use nextClear and nextSet check interval
			clearIndex = chrBitSet.nextClearBit(fromIndex);
		}

	}
}
