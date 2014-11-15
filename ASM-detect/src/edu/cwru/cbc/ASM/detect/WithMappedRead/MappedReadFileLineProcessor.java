package edu.cwru.cbc.ASM.detect.WithMappedRead;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/21/14.
 * Implemented LineProcessor for reading MappedRead File
 * Mapped reads start pos is 1-based, end pos is 0-based.
 */
public class MappedReadFileLineProcessor implements LineProcessor<List<MappedRead>> {
	private List<MappedRead> mappedReadList = new ArrayList<>();

	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr") || line.startsWith("ref")) {
			return true;
		} else if (line.equals("")) {
			return false;
		} else {
			String[] items = line.split("\t");
			if (items[1].length() != 1) {
				throw new RuntimeException("invalid strand!");
			}
			mappedReadList.add(
					new MappedRead(items[0], items[1].charAt(0), Integer.parseInt(items[2]), Integer.parseInt(items[3]),
								   items[4], items[5]));
			return true;
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return mappedReadList;
	}
}
