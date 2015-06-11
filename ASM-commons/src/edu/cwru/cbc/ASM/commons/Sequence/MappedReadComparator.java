package edu.cwru.cbc.ASM.commons.Sequence;

import java.util.Comparator;

//TODO add unit test
public class MappedReadComparator implements Comparator<MappedRead> {

	@Override
	public int compare(MappedRead r1, MappedRead r2) {
		// rank smaller start, longer end in front
		if (r1.getStart() < r2.getStart()) {
			return -1;
		} else if (r1.getStart() > r2.getStart()) {
			return 1;
		} else if (r1.getEnd() < r2.getEnd()) {
			return -1;
		} else if (r1.getEnd() > r2.getEnd()) {
			return 1;
		} else {
			return 0;
		}
	}
}
