package edu.cwru.cbc.ASM.detect.Utils;

import edu.cwru.cbc.ASM.detect.DataType.MappedRead;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/21/14.
 * Implemented LineProcessor for reading MappedRead File
 */
public class MappedReadFileLineProcessor implements LineProcessor<List<MappedRead>> {
	private List<MappedRead> mappedReadList = new ArrayList<>();

	@Override
	public boolean processLine(String line) throws IOException {
        if (line.startsWith("chr") || line.startsWith("ref")) {
            return true;
		}else if(line.equals("")){
			return false;
		}else {
			List<String> itemList = Lists.newArrayList(Splitter.on('\t').split(line));
			mappedReadList.add(new MappedRead(itemList.get(0), itemList.get(1), Long.parseLong(itemList.get(2)), Long.parseLong(itemList.get(3)), itemList.get(4),
					Long.parseLong(itemList.get(5))));
			return true;
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return mappedReadList;
	}
}
