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
public class GroupedReadsLineProcessor implements LineProcessor<Pair<String, List<List<MappedRead>>>> {
	private static final Splitter tabSplitter = Splitter.on("\t");
	private String ref;
	private List<List<MappedRead>> groups = new ArrayList<>();
	private List<MappedRead> group = new ArrayList<>();

	@Override
	public boolean processLine(@Nonnull String line) throws IOException {
		List<String> itemList = tabSplitter.splitToList(line);
		if (line.startsWith("ref:")) {
			ref = itemList.get(1);
		} else if (line.equals("")) {
			groups.add(group);
			group = new ArrayList<>();
		} else {
			// TODO add input data validation
			group.add(new MappedRead(itemList.get(0), itemList.get(1).charAt(0), Integer.parseInt(itemList.get(2)),
					itemList.get(1).charAt(0) == '+' ?
							itemList.get(4).replace(".", "") :
							MappedRead.getComplementarySequence(itemList.get(4).replace(".", "")),
					itemList.get(5)));
		}
		return true;
	}

	@Override
	public Pair<String, List<List<MappedRead>>> getResult() {
		return new ImmutablePair<>(ref, groups);
	}

}
