package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Created by kehu on 11/4/16.
 * SAM/BAM file reader using htsjdk API.
 */
public class MappedReadReader {
	public static List<MappedRead> readMappedReads(File inputSamOrBamFile) throws IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
		List<MappedRead> list = readMappedReads(reader.iterator());
		reader.close();
		return list;
	}

	public static List<MappedRead> readMappedReads(File inputSamOrBamFile, String chr) throws IOException {
		return readMappedReads(inputSamOrBamFile, chr, 0, 0);
	}

	public static List<MappedRead> readMappedReads(File inputSamOrBamFile, String chr, int start, int end) throws
			IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
		SAMRecordIterator iterator = reader.query(chr, start, end, false);
		List<MappedRead> list = readMappedReads(iterator);
		iterator.close();
		reader.close();
		return list;
	}

	private static List<MappedRead> readMappedReads(SAMRecordIterator iterator) {
		LinkedHashMap<String, MappedRead> mappedReadLinkedHashMap = new LinkedHashMap<>();
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			MappedRead mappedRead = new MappedRead(samRecord.getReferenceName(),
					samRecord.getReadNegativeStrandFlag() ? '-' : '+', samRecord.getStart(), samRecord.getReadString(),
					samRecord.getReadName(), samRecord.getFirstOfPairFlag());
			if (samRecord.getReadPairedFlag()) {
				if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
					MappedRead existRead = mappedReadLinkedHashMap.get(mappedRead.getId());
					if (isInOrder(existRead, mappedRead)) {
						// existRead is first
						existRead.updateSequence(combinePERead(existRead.getSequence(), mappedRead.getSequence(),
								mappedRead.getEnd() - existRead.getStart() + 1));
					} else {
						// mappedRead is first
						existRead.updateSequence(combinePERead(mappedRead.getSequence(), existRead.getSequence(),
								existRead.getEnd() - mappedRead.getStart() + 1));
						existRead.setStart(mappedRead.getStart());
					}
				} else {
					mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
				}
			} else {
				if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
					throw new RuntimeException("found duplicate mapped read!");
				} else {
					mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
				}
			}
		}
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

	private static boolean isInOrder(MappedRead existRead, MappedRead newRead) {
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

	private static String combinePERead(String p1, String p2, int length) {
		if (length < 0) {
			throw new RuntimeException("invalid new length:" + length);
		}
		StringBuilder sb = new StringBuilder(length);
		int offset = length - p1.length() - p2.length();
		if (offset < 0) {
			sb.append(p1).append(p2.substring(-1 * offset));
		} else {
			sb.append(p1).append(StringUtils.repeat('-', offset)).append(p2);
		}
		return sb.toString();
	}
}
