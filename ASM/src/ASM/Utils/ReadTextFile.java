package ASM.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/14/14.
 */
public class ReadTextFile {
	public static <T> List<T> ReadTextFileWithDelimiter(Constructor<T> constructor, String fileName, String delimiter, boolean hasColumnName) throws IOException {
		List<T> resultList = new ArrayList<>();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		String line;
		String[] items;
		if (hasColumnName) {
			bufferedReader.readLine();
		}
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split(delimiter);
			Class<?>[] constructorTypes = constructor.getParameterTypes();
			Object[] objects = new Object[constructorTypes.length];
			// assume items have same order with constructor parameters
			for (int i = 0; i < constructorTypes.length; i++) {
				try {
					if (constructorTypes[i].isPrimitive()) {
						objects[i] = ReflectionUtils.parseObjectFromString(items[i], constructorTypes[i]);
					}else {
						objects[i] = ReflectionUtils.parseObjectFromString(items[i], constructorTypes[i]);
					}
				} catch (NoSuchMethodException | InvocationTargetException | InstantiationException | IllegalAccessException
						e) {
					e.printStackTrace();
				}
			}
			try {
				resultList.add(constructor.newInstance(objects));
			} catch (InstantiationException | IllegalAccessException | InvocationTargetException e) {
				e.printStackTrace();
			}
		}

		return resultList;
	}
}
