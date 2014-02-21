package ASM.Utils;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.BitSet;
import java.util.List;

/**
 * Created by ke on 2/21/14.
 * Implement LineProcessor for generating interval summary
 */
public class IntervalChekingLineProcessor implements LineProcessor<BitSet> {
	private BitSet chrBitSet;
//	private int[] chrCounter;
	private int counter=0;

	public IntervalChekingLineProcessor(int chrBitSize) {
		this.chrBitSet = new BitSet(chrBitSize);
//		this.chrCounter = new int[chrBitSize];
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr")) {
			return true;
		}else if (line.equals("")){
			return false;
		}else{
			List<String> itemList = Lists.newArrayList(Splitter.on('\t').split(line));
			chrBitSet.set(Integer.parseInt(itemList.get(2)), Integer.parseInt(itemList.get(3)));
			counter++;
			if (counter % 1000000 == 0) {
				System.out.println(counter);
			}
			return true;
		}
	}

	@Override
	public BitSet getResult() {
		return null;
	}
}
