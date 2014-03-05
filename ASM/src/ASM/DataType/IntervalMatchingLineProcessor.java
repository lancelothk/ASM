package ASM.DataType;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.List;

/**
 * Created by ke on 3/5/14.
 * Implement interval and read matching
 */
public class IntervalMatchingLineProcessor implements LineProcessor<List<GenomicInterval>> {
	private List<GenomicInterval> intervalList;
	private int counter = 0;

	public IntervalMatchingLineProcessor(List<GenomicInterval> intervalList) {
		System.out.println("start interval matching:");
		this.intervalList = intervalList;
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr")) {
			return true;
		} else if (line.equals("")) {
			return false;
		} else {
			List<String> itemList = Lists.newArrayList(Splitter.on('\t').split(line));
			int start = Integer.parseInt(itemList.get(2));
			int end = Integer.parseInt(itemList.get(3));
			for (GenomicInterval interval : intervalList) {
				if (start <= interval.getStart() && interval.getEnd() <= end) {
					interval.incrementReadCount();
				}
			}
			counter++;
			if (counter % 1000000 == 0) {
				System.out.println(counter);
			}
			return true;
		}
	}

	@Override
	public List<GenomicInterval> getResult() {
		return this.intervalList;
	}
}
