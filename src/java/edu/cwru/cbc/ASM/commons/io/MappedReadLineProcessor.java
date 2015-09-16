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
	protected LinkedHashMap<String, MappedRead> mappedReadLinkedHashMap = new LinkedHashMap<>();
	protected Predicate<MappedRead> criteria;
	private boolean isPairEnd;

	public MappedReadLineProcessor() {
		this(false, mr -> true);

	}

	public MappedReadLineProcessor(boolean isPairEnd, Predicate<MappedRead> criteria) {
		this.criteria = criteria;
		this.isPairEnd = isPairEnd;
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr\t") || line.startsWith("ref") || line.startsWith("assembly")) {
			return true;
		} else {
			List<String> itemList = tabSplitter.splitToList(line);
			if (isPairEnd) {
				if (!itemList.get(1).equals("+") && !itemList.get(1).equals("-")) {
					throw new RuntimeException("invalid strand! in line:\t" + line);
				}
				if (!IUPACCode.validateNucleotideCode(itemList.get(4)) || !IUPACCode.validateNucleotideCode(
						itemList.get(5))) {
					throw new RuntimeException("invalid character in sequence!\t" + line);
				}
				addMappedReadToMap(line, new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
						Integer.parseInt(itemList.get(2)),
						itemList.get(4),
						itemList.get(6) + "_1"));
				addMappedReadToMap(line, new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
						Integer.parseInt(itemList.get(3)) - itemList.get(5).length(),
						itemList.get(5),
						itemList.get(6) + "_2"));
			} else {
				if (itemList.size() == 6) {
					if (!itemList.get(1).equals("+") && !itemList.get(1).equals("-")) {
						throw new RuntimeException("invalid strand! in line:\t" + line);
					}
					if (!IUPACCode.validateNucleotideCode(itemList.get(4))) {
						throw new RuntimeException("invalid character in sequence!\t" + line);
					}
					// h1, i90
					addMappedReadToMap(line, new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
							Integer.parseInt(itemList.get(2)),
							itemList.get(4),
							itemList.get(5)));
				} else if (itemList.size() == 10) {
					if (!itemList.get(4).equals("+") && !itemList.get(4).equals("-")) {
						throw new RuntimeException("invalid strand! in line:\t" + line);
					}
					if (!IUPACCode.validateNucleotideCode(itemList.get(8))) {
						throw new RuntimeException("invalid character in sequence!\t" + line);
					}
					// ff, h9_laurent
					addMappedReadToMap(line, new MappedRead(itemList.get(3), itemList.get(4).charAt(0),
							Integer.parseInt(itemList.get(5)),
							itemList.get(8),
							itemList.get(6)));
				} else {
					throw new RuntimeException("incompatible format detected! column number is " + itemList.size());
				}
			}
			return true;
		}
	}

	private void addMappedReadToMap(String line, MappedRead mappedRead) {
		if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
			throw new RuntimeException("found duplicate mapped read! in line:\t" + line);
		} else if (criteria.test(mappedRead)) {
			mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

}
