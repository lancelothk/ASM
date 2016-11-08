package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
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
	public static void readMappedReads(File inputSamOrBamFile, HashIntObjMap<RefCpG> refMap,
	                                   MappedReadHandler handler, String chr, int start, int end) throws IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
		SAMRecordIterator iterator = reader.query(chr, start, end, false);
		readMappedReads(iterator, refMap, handler);
		iterator.close();
		reader.close();
	}

	private static void readMappedReads(SAMRecordIterator iterator, HashIntObjMap<RefCpG> refMap,
	                                    MappedReadHandler handler) throws IOException {
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			MappedRead mappedRead = new MappedRead(samRecord.getReferenceName(),
					samRecord.getReadNegativeStrandFlag() ? '-' : '+', samRecord.getStart(),
					samRecord.getReadString(), samRecord.getReadName());
			mappedRead.generateCpGsInRead(refMap);
			handler.processMappedRead(mappedRead);
		}
	}

	public static List<MappedRead> readMappedReads(File inputSamOrBamFile, HashIntObjMap<RefCpG> refMap,
	                                               boolean isPairedEnd, String chr) throws IOException {
		return readMappedReads(inputSamOrBamFile, refMap, isPairedEnd, chr, 0, 0);
	}

	public static List<MappedRead> readMappedReads(File inputSamOrBamFile, HashIntObjMap<RefCpG> refMap,
	                                               boolean isPairedEnd, String chr, int start,
	                                               int end) throws IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
		SAMRecordIterator iterator = reader.query(chr, start, end, false);
		List<MappedRead> list = readMappedReads(iterator, refMap, isPairedEnd);
		iterator.close();
		reader.close();
		return list;
	}

	private static List<MappedRead> readMappedReads(SAMRecordIterator iterator, HashIntObjMap<RefCpG> refMap,
	                                                boolean isPairedEnd) {
		LinkedHashMap<String, MappedRead> mappedReadLinkedHashMap = new LinkedHashMap<>();
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			if (isPairedEnd) {
				MappedRead mappedRead = new MappedRead(samRecord.getReferenceName(),
						samRecord.getReadNegativeStrandFlag() ? '-' : '+', samRecord.getStart(),
						samRecord.getReadString(),
						samRecord.getReadName(), samRecord.getFirstOfPairFlag());
				if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
					MappedRead existRead = mappedReadLinkedHashMap.get(mappedRead.getId());
					if (existRead.getStrand() != mappedRead.getStrand()) {
						System.err.println("invalid paired end data!(from different strand):" + existRead.getId());
						break;
					}
					if (existRead.isFirstInPair() == mappedRead.isFirstInPair()) {
						throw new RuntimeException(
								"invalid paired end data!(have same flag):" + existRead.getId() + "-" + mappedRead.getId());
					}
					if (isInOrder(existRead)) {
						// existRead is first
						existRead.updateSequence(combinePERead(existRead.getSequence(), mappedRead.getSequence(),
								mappedRead.getEnd() - existRead.getStart() + 1));
					} else {
						// mappedRead is first
						existRead.updateSequence(combinePERead(mappedRead.getSequence(), existRead.getSequence(),
								existRead.getEnd() - mappedRead.getStart() + 1));
						existRead.setStart(mappedRead.getStart());
					}
					existRead.generateCpGsInRead(refMap);
				} else {
					mappedRead.generateCpGsInRead(refMap);
					mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
				}
			} else {
				MappedRead mappedRead = new MappedRead(samRecord.getReferenceName(),
						samRecord.getReadNegativeStrandFlag() ? '-' : '+', samRecord.getStart(),
						samRecord.getReadString(), samRecord.getReadName());
				if (mappedReadLinkedHashMap.containsKey(mappedRead.getId())) {
					throw new RuntimeException("found duplicate mapped read!");
				} else {
					mappedRead.generateCpGsInRead(refMap);
					mappedReadLinkedHashMap.put(mappedRead.getId(), mappedRead);
				}
			}
		}
		return new ArrayList<>(mappedReadLinkedHashMap.values());
	}

	private static boolean isInOrder(MappedRead existRead) {
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
