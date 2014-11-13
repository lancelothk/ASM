package edu.cwru.cbc.ASM.align;

import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by ke on 2/21/14.
 * Test Read.toString()
 */
public class ReadTest {
	@Test
	public void testToString() throws Exception {
		MappedRead read = new MappedRead("6", '+', 10, 20, "1234567890", "test");
		assertEquals(read.toString(5), "6\t+\t10\t20\t.....1234567890\ttest");
		System.out.println(read.toString(5));
	}
}
