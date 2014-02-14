package ASM.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/14/14.
 */
public class ReadTextFile {
	public <T> List<T> ReadTextFileWithDelimiter(Class<T> clazz, String fileName, String delimiter, boolean hasColumnName) throws IOException {
		List<T> resultList = new ArrayList<>();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		String line;
		String[] items;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split(delimiter);

		}
		return resultList;
	}
}
