package edu.cwru.cbc.ASM.commons.ReflectionUtils;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

/**
 * Created by kehu on 12/9/14.
 * Utils for using reflection
 */
public class ReflectionUtils {
    public static void setFinalStatic(Field field, Object newValue) throws Exception {
        field.setAccessible(true);

        Field modifiersField = Field.class.getDeclaredField("modifiers");
        modifiersField.setAccessible(true);
        modifiersField.setInt(field, field.getModifiers() & ~Modifier.FINAL);

        field.set(null, newValue);
    }

    public static void invokePrivateMethod(Method method, Object obj,
                                           Object... args) throws InvocationTargetException, IllegalAccessException {
        method.setAccessible(true);
        method.invoke(obj, args);
    }
}
