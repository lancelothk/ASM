package edu.cwru.cbc.ASM.commons.io;

import com.google.common.base.Splitter;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import javax.annotation.Nonnull;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 7/6/16.
 * Class to read grouped reads file.
 */
public class GroupedReadsLineProcessor implements LineProcessor<Pair<String, Pair<List<MappedRead>, List<MappedRead>>>> {
	private static final Splitter tabSplitter = Splitter.on("\t");
	private String ref;
	private int group = 0;
	private List<MappedRead> group1 = new ArrayList<>();
	private List<MappedRead> group2 = new ArrayList<>();

	@Override
	public boolean processLine(@Nonnull String line) throws IOException {
		List<String> itemList = tabSplitter.splitToList(line);
		if (line.startsWith("ref:")) {
			ref = itemList.get(1);
		} else if (line.equals("")) {
			group++;
		} else {
			// TODO add input data validation
			switch (group) {
				case 0:
					group1.add(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
							Integer.parseInt(itemList.get(2)),
							itemList.get(1).charAt(0) == '+' ? itemList.get(4)
									.replace(".", "") : MappedRead.getComplementarySequence(
									itemList.get(4).replace(".", "")), itemList.get(5)));
					break;
				case 1:
					group2.add(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
							Integer.parseInt(itemList.get(2)),
							itemList.get(1).charAt(0) == '+' ? itemList.get(4)
									.replace(".", "") : MappedRead.getComplementarySequence(
									itemList.get(4).replace(".", "")), itemList.get(5)));
					break;
				default:
					throw new RuntimeException("more than 2 groups!");
			}
		}
		return true;
	}

	@Override
	public Pair<String, Pair<List<MappedRead>, List<MappedRead>>> getResult() {
		return new ImmutablePair<>(ref, new ImmutablePair<>(group1, group2));
	}

}
