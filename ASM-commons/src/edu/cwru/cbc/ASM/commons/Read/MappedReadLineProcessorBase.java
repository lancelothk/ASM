package edu.cwru.cbc.ASM.commons.Read;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.CpG.CpG;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.MethylStatus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by lancelothk on 6/8/15.
 * abstract base class for MappedReadLineProcessor.
 * Shared by both MappedReadLineProcessor and MappedReadLineProcessorWithSummary
 * Mapped reads start pos is 1-based, end pos is 0-based.
 */
public abstract class MappedReadLineProcessorBase<T> implements LineProcessor<T> {
	protected List<MappedRead> mappedReadList = new ArrayList<>();
	protected Map<Integer, RefCpG> refMap;

	// don't filter read with min_read_cpg threshold
	public MappedReadLineProcessorBase(List<RefCpG> refCpGList) {
		this.refMap = refCpGList.stream().collect(Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
	}


	@Override
	public boolean processLine(String line) throws IOException {
		try {
			if (line.startsWith("chr") || line.startsWith("ref") || line.startsWith("assembly")) {
				return true;
			} else if (line.equals("")) {
				return false;
			} else {
				updateRefCpG(processRead(line));
				return true;
			}
		} catch (Exception e) {
			throw new RuntimeException("Problem line: " + line + "\n", e);
		}
	}

	protected MappedRead processRead(String line) {
		String[] items = line.split("\t");
		if (!items[1].equals("+") && !items[1].equals("-")) {
			throw new RuntimeException("invalid strand!");
		}
		int start = Integer.parseInt(items[2]);// mapped read is 1bp right than UCSC ref(1 based).
		int end = Integer.parseInt(items[3]);

		MappedRead mappedRead;
		if (items.length == 6) {
			// for h1/i90 dataset
			mappedRead = new MappedRead(items[0], items[1].charAt(0), start, end, items[4], items[5]);
		} else if (items.length == 7) {
			// for h9 and other dataset
			mappedRead = new MappedRead(items[0], items[1].charAt(0), start, end, items[4], items[6]);
		} else {
			throw new RuntimeException("columns is not correct for mapped read format");
		}
		return mappedRead;
	}

	protected void updateRefCpG(MappedRead mappedRead) {
		for (int i = mappedRead.getStart(); i < mappedRead.getEnd(); i++) {
			if (refMap.containsKey(i)) {
				CpG cpg = new CpG(mappedRead, refMap.get(i), mappedRead.getMethylStatus(i));
				// ignore unknown methyl CpG
				if (cpg.getMethylStatus() != MethylStatus.N) {
					mappedRead.addCpG(cpg);
					refMap.get(i).addCpG(cpg);
				}
				i++;
			}
		}
		mappedReadList.add(mappedRead);
	}
}
