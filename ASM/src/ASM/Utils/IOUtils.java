package ASM.Utils;

import ASM.DataType.MappedRead;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/19/14.
 */
public class IOUtils {
	public static List<MappedRead> readMappedRead(String fileName) throws IOException {
		List<MappedRead> mappedReadList = new ArrayList<>();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		String line;
		String[] items;
		while ((line = bufferedReader.readLine()) != null){
			items = line.split("\t");
			MappedRead mappedRead = new MappedRead(items[0], items[1], Long.parseLong(items[2]), Long.parseLong(items[3]), items[4],
					Long.parseLong(items[5]));
			mappedReadList.add(mappedRead);
		}
		bufferedReader.close();
		return mappedReadList;
	}
}
