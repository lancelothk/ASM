package ASM.Utils;

import ASM.DataType.MappedRead;
import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/21/14.
 */
public class MappedReadFileLineProcessor implements LineProcessor<List<MappedRead>> {
	List<MappedRead> mappedReadList = new ArrayList<>();
	@Override
	public boolean processLine(String line) throws IOException {
		if (line.startsWith("chr")){
			return true;
		}else if(!line.equals("")){
			String[] items = line.split("\t");
			MappedRead mappedRead = new MappedRead(items[0], items[1], Long.parseLong(items[2]), Long.parseLong(items[3]), items[4],
					Long.parseLong(items[5]));
			mappedReadList.add(mappedRead);
			return true;
		}else {
			return false;
		}
	}

	@Override
	public List<MappedRead> getResult() {
		return mappedReadList;
	}
}
