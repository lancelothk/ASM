package edu.cwru.cbc.ASM.commons.io;

import com.google.common.base.Splitter;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.sequence.IUPACCode;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Predicate;

/**
 * Created by lancelothk on 6/8/15.
 * abstract base class for MappedReadLineProcessor.
 * Shared by both MappedReadLineProcessor and MappedReadLineProcessorWithSummary
 * Mapped reads start pos is 0-based, end pos is 0-based.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
	private static final Splitter tabSplitter = Splitter.on("\t");
	protected LinkedHashMap<Integer, MappedRead> mappedReadLinkedHashMap = new LinkedHashMap<>();
	protected Predicate<MappedRead> criteria;

	public MappedReadLineProcessor() {
		criteria = mr -> true;
	}

	public MappedReadLineProcessor(Predicate<MappedRead> criteria) {
		this.criteria = criteria;
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr") || line.startsWith("ref") || line.startsWith("assembly")) {
			return true;
		} else {
			MappedRead mr = processRead(line);
			if (mappedReadLinkedHashMap.containsKey(mr.getId())) {
				throw new RuntimeException("found duplicate mapped read! in line:\t" + line);
			} else if (criteria.test(mr)) {
				mappedReadLinkedHashMap.put(mr.getId(), mr);
			}
			return true;
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

	private MappedRead processRead(String line) {
		List<String> itemList = tabSplitter.splitToList(line);
		if (itemList.size() != 6 & itemList.size() != 7) {
			throw new RuntimeException("columns is not correct for mapped read format in line:\t" + line);
		}
		if (!itemList.get(1).equals("+") && !itemList.get(1).equals("-")) {
			throw new RuntimeException("invalid strand! in line:\t" + line);
		}
		int start = Integer.parseInt(itemList.get(2));// mapped read is 0 based start.

		if (!IUPACCode.validateNucleotideCode(itemList.get(4))) {
			throw new RuntimeException("invalid character in sequence!\t" + line);
		}
		MappedRead mappedRead;
		if (itemList.size() == 6) {
			mappedRead = new MappedRead(itemList.get(0), itemList.get(1).charAt(0), start, itemList.get(4),
					Integer.parseInt(itemList.get(5)));
		} else {
			mappedRead = new MappedRead(itemList.get(0), itemList.get(1).charAt(0), start, itemList.get(4),
					Integer.parseInt(itemList.get(6)));
		}
		return mappedRead;
	}

}
