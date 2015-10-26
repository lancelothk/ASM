package edu.cwru.cbc.ASM.commons.io;

import com.google.common.base.Splitter;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.MappedReadFileFormat;
import edu.cwru.cbc.ASM.commons.sequence.IUPACCode;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.apache.commons.lang3.StringUtils;

import javax.annotation.Nonnull;
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
	public boolean processLine(@Nonnull String line) throws IOException {
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
				processSAMRead(itemList, isPairEnd);
				return true;
			}
		} catch (Exception e) {
			System.err.println(line);
			throw e;
		}
	}

	private void processSAMRead(List<String> itemList, boolean isPairEnd) {
		char strand;
		if (isPairEnd) {
			if (Integer.parseInt(itemList.get(1)) == 99 || Integer.parseInt(itemList.get(1)) == 147) {
				strand = '+';
			} else if (Integer.parseInt(itemList.get(1)) == 83 || Integer.parseInt(itemList.get(1)) == 163) {
				strand = '-';
			} else {
				throw new RuntimeException("invalid paired end bismark sam flag!");
			}
			validateRead(itemList.get(9));
			// sam format position is 1-based
			addPairedMappedReadToMap(
					new MappedRead(itemList.get(2), strand, Integer.parseInt(itemList.get(3)) - 1,
							strand == '-' ? MappedRead.getComplementarySequence(itemList.get(9)) : itemList.get(9),
							itemList.get(0),
							Integer.parseInt(itemList.get(1)) == 99 || Integer.parseInt(itemList.get(1)) == 83));
		} else {
			if (Integer.parseInt(itemList.get(1)) == 0) {
				strand = '+';
			} else if (Integer.parseInt(itemList.get(1)) == 16) {
				strand = '-';
			} else {
				throw new RuntimeException("invalid single end bismark sam flag!");
			}
			validateRead(itemList.get(9));
			// sam format position is 1-based
			addMappedReadToMap(
					new MappedRead(itemList.get(2), strand, Integer.parseInt(itemList.get(3)) - 1,
							strand == '-' ? MappedRead.getComplementarySequence(itemList.get(9)) : itemList.get(9),
							itemList.get(0)));
		}
	}

	private void addPairedMappedReadToMap(MappedRead newRead) {
		if (mappedReadLinkedHashMap.containsKey(newRead.getId())) {
			MappedRead existRead = mappedReadLinkedHashMap.get(newRead.getId());
			if (isInOrder(existRead, newRead)) {
				// existRead is first
				existRead.updateSequence(combinePERead(existRead.getSequence(), newRead.getSequence(),
						newRead.getEnd() - existRead.getStart() + 1));
			} else {
				// mappedRead is first
				existRead.updateSequence(combinePERead(newRead.getSequence(), existRead.getSequence(),
						existRead.getEnd() - newRead.getStart() + 1));
				existRead.setStart(newRead.getStart());
			}
			if (!criteria.test(existRead)) {
				mappedReadLinkedHashMap.remove(existRead.getId());
			}
		} else {
			mappedReadLinkedHashMap.put(newRead.getId(), newRead);
		}
	}

	private boolean isInOrder(MappedRead existRead, MappedRead newRead) {
		if (existRead.getStrand() != newRead.getStrand()) {
			throw new RuntimeException("invalid paired end data!(from different strand)");
		}
		if (existRead.isFirstInPair() == newRead.isFirstInPair()) {
			throw new RuntimeException("invalid paired end data!(have same flag)");
		}
		if (existRead.getStrand() == '+') {
			return existRead.isFirstInPair();
		} else if (existRead.getStrand() == '-') {
			return !existRead.isFirstInPair();
		} else {
			throw new RuntimeException("invalid strand!");
		}
	}

	private String combinePERead(String p1, String p2, int length) {
		if (length < 0) {
			throw new RuntimeException("invalid new length:" + length);
		}
		@SuppressWarnings("StringBufferReplaceableByString")
		StringBuilder sb = new StringBuilder(length);
		int offset = length - p1.length() - p2.length();
		if (offset < 0) {
			sb.append(p1).append(p2.substring(-1 * offset));
		} else {
			sb.append(p1).append(StringUtils.repeat('-', offset)).append(p2);
		}
		return sb.toString();
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
							Integer.parseInt(itemList.get(2)), combinePERead(itemList.get(4), itemList.get(5),
							Integer.parseInt(itemList.get(3)) - Integer.parseInt(itemList.get(2))), itemList.get(6)));
				} else {
					if (itemList.size() == 6) {
						validateStrand(itemList.get(1));
						validateRead(itemList.get(4));
						// h1, i90
						addMappedReadToMap(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
								Integer.parseInt(itemList.get(2)), itemList.get(4), itemList.get(5)));
					} else if (itemList.size() == 7) {
						validateStrand(itemList.get(1));
						validateRead(itemList.get(4));
						// h9
						addMappedReadToMap(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
								Integer.parseInt(itemList.get(2)), itemList.get(4), itemList.get(6)));
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
