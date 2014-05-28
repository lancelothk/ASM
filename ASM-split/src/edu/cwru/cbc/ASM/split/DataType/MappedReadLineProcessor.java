package edu.cwru.cbc.ASM.split.DataType;

import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by lancelothk on 5/27/14.
 * LineProcessor for process mapped read.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
	private Map<Integer, CpGSite> refMap;
	private List<MappedRead> mappedReadList;

	public MappedReadLineProcessor(Map<Integer, CpGSite> refMap) {
		this.refMap = refMap;
		this.mappedReadList = new ArrayList<>();
	}

	@Override
	public boolean processLine(String s) throws IOException {
		String[] items = s.split("\t");
		if (items.length != 6){
			throw new RuntimeException("invalid mapped read format:" + s);
		}
		MappedRead mappedRead = new MappedRead(items[0], items[1], items[5]);
		//to suit prev program
		mappedRead.setSequence(items[4]);
		mappedRead.setStart(Integer.parseInt(items[2]));
		mappedRead.setEnd(Integer.parseInt(items[3]));

		int start = Integer.parseInt(items[2]) - 1;// mapped read is 1bp right than UCSC ref
		int end = Integer.parseInt(items[3]) - 1;
		for (int i = start; i < end; i++) {
			if (refMap.containsKey(i)){
				CpG cpg = new CpG(mappedRead, refMap.get(i));
				cpg.setMethylated(extractMethylStatus(items[1], items[4].substring(i - start, i - start + 2)));
				mappedRead.addCpG(cpg);
				refMap.get(i).addCpG(cpg);
				i++;
			}
		}
		mappedReadList.add(mappedRead);
		return true;
	}

	private boolean extractMethylStatus(String strand, String sequence) {
		if (sequence.length() != 2){
			throw new RuntimeException("invalid input cpg sequence:\t" + sequence);
		}
		// only consider 'C' position in CpG. Ignore the char in 'G' position
		switch (strand) {
			case "+": {
				// CN is methylated, TN is non-methylated
				char bp = sequence.charAt(0);
				if (bp == 'C') {
					return true;
				} else if (bp == 'T') {
					return false;
				} else {
					// TODO should ignore non C/T case in CpG list
					return false;
				}
			}
			case "-": {
				// NC is methylated, NT is non-methylated
				char bp = sequence.charAt(1);
				if (bp == 'C') {  // for complementary bp, 'G'
					return true;
				} else if (bp == 'T') { // for complementary bp, 'A'
					return false;
				} else {
					// TODO should ignore non C/T case in CpG list
					return false;
				}
			}
			default:
				throw new RuntimeException("illegal strand!");
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return this.mappedReadList;
	}
}
