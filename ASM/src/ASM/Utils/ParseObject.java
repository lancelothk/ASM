package ASM.Utils;

import java.lang.reflect.Method;

/**
 * Created by lancelothk on 2/15/14.
 */
public interface ParseObject<T> {
	public T parse(Method constructor);
}
