package ASM.Utils;

import ASM.DataType.ChrCoverageSummary;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.List;

/**
 * Created by ke on 2/21/14.
 * Implement LineProcessor for generating interval summary
 */
public class IntervalChekingLineProcessor implements LineProcessor<ChrCoverageSummary> {
	private ChrCoverageSummary chrCoverageSummary;
	private int counter=0;

	public IntervalChekingLineProcessor(int chrBitSize, String referenceFileName) throws IOException {
		this.chrCoverageSummary = new ChrCoverageSummary(chrBitSize, referenceFileName);
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr")) {
			return true;
		}else if (line.equals("")){
			return false;
		}else{
			List<String> itemList = Lists.newArrayList(Splitter.on('\t').split(line));
			int start = Integer.parseInt(itemList.get(2));
			int end = Integer.parseInt(itemList.get(3));
			this.chrCoverageSummary.addCoverage(start, end);
			counter++;
			if (counter % 1000000 == 0) {
				System.out.println(counter);
			}
			return true;
		}
	}

	@Override
	public ChrCoverageSummary getResult() {
		return this.chrCoverageSummary;
	}
}
