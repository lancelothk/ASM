package edu.cwru.cbc.ASM.commons.io;

import com.google.common.base.Splitter;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.MappedReadFileFormat;
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
	private MappedReadFileFormat inputFormat;

	public MappedReadLineProcessor() {
		this(false, mr -> true, MappedReadFileFormat.MAPPEDREAD);
	}

	public MappedReadLineProcessor(Predicate<MappedRead> criteria) {
		this(false, criteria, MappedReadFileFormat.MAPPEDREAD);
	}

	public MappedReadLineProcessor(boolean isPairEnd, Predicate<MappedRead> criteria) {
		this(isPairEnd, criteria, MappedReadFileFormat.MAPPEDREAD);
	}

	public MappedReadLineProcessor(boolean isPairEnd, Predicate<MappedRead> criteria,
	                               MappedReadFileFormat inputFormat) {
		this.criteria = criteria;
		this.isPairEnd = isPairEnd;
		this.inputFormat = inputFormat;
	}

	@Override
	public boolean processLine(String line) throws IOException {
		if (inputFormat == MappedReadFileFormat.MAPPEDREAD) {
			return processMappedReadFormat(line);
		} else if (inputFormat == MappedReadFileFormat.SAM) {
			return processSAMFormat(line);
		} else {
			throw new RuntimeException("unknown input format!");
		}
	}

	private boolean processSAMFormat(String line) {
		try {
			if (line.startsWith("@")) {
				return true;
			} else {
				List<String> itemList = tabSplitter.splitToList(line);
				if (isPairEnd) {
					//TODO
				} else {
					char strand;
					if ((Integer.parseInt(itemList.get(1)) & 16) == 16) {
						strand = '-';
					} else {
						strand = '+';
					}
					validateRead(itemList.get(9));
					// sam format position is 1-based
					addMappedReadToMap(
							new MappedRead(itemList.get(2), strand, Integer.parseInt(itemList.get(3)) - 1,
									itemList.get(9),
									itemList.get(0)));
				}
				return true;
			}
		} catch (Exception e) {
			System.err.println(line);
			throw e;
		}
	}

	private boolean processMappedReadFormat(String line) {
		try {
			if (line.startsWith("chr\t") || line.startsWith("ref") || line.startsWith("assembly")) {
				return true;
			} else {
				List<String> itemList = tabSplitter.splitToList(line);
				if (isPairEnd) {
					validateStrand(itemList.get(1));
					validateRead(itemList.get(4));
					validateRead(itemList.get(5));
					addMappedReadToMap(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
							Integer.parseInt(itemList.get(2)), itemList.get(4), itemList.get(6) + "_1"));
					addMappedReadToMap(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
							Integer.parseInt(itemList.get(3)) - itemList.get(5).length(), itemList.get(5),
							itemList.get(6) + "_2"));
				} else {
					if (itemList.size() == 6) {
						validateStrand(itemList.get(1));
						validateRead(itemList.get(4));
						// h1, i90
						addMappedReadToMap(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
								Integer.parseInt(itemList.get(2)), itemList.get(4), itemList.get(5)));
					} else if (itemList.size() == 10) {
						validateStrand(itemList.get(4));
						validateRead(itemList.get(8));
						// ff, h9_laurent
						addMappedReadToMap(new MappedRead(itemList.get(3), itemList.get(4).charAt(0),
								Integer.parseInt(itemList.get(5)), itemList.get(8), itemList.get(6)));
					} else {
						throw new RuntimeException("incompatible format detected! column number is " + itemList.size());
					}
				}
				return true;
			}
		} catch (Exception e) {
			System.err.println(line);
			throw e;
		}
	}

	private void validateRead(String sequence) {
		if (!IUPACCode.validateNucleotideCode(sequence)) {
			throw new RuntimeException("invalid character in sequence!");
		}
	}

	private void validateStrand(String strand) {
		if (!strand.equals("+") && !strand.equals("-")) {
			throw new RuntimeException("invalid strand!");
		}
	}

	private void addMappedReadToMap(MappedRead mappedRead) {
		if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
			throw new RuntimeException("found duplicate mapped read!");
		} else if (criteria.test(mappedRead)) {
			mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

}
