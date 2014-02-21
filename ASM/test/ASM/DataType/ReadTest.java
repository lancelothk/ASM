package ASM.DataType;

import junit.framework.Assert;
import org.junit.Test;

/**
 * Created by ke on 2/21/14.
 */
public class ReadTest {
	@Test
	public void testToString() throws Exception {
		Read read = new Read("6", '+',10,20,"1234567890","test");
		Assert.assertEquals(read.toString(5),"6\t+\t10\t20\t.....1234567890\ttest");
		System.out.println(read.toString(5));
	}
}
