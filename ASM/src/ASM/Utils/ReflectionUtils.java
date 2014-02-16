package ASM.Utils;

import java.lang.reflect.InvocationTargetException;

/**
 * Created by lancelothk on 2/16/14.
 */
public class ReflectionUtils {
	public static <T> T parseObjectFromString(String s, Class<T> clazz) throws NoSuchMethodException, IllegalAccessException, InvocationTargetException, InstantiationException {
		return clazz.getConstructor(new Class[]{String.class}).newInstance(s);
	}
}
