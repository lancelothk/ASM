package edu.cwru.cbc.ASM.commons.Sequence;

import edu.cwru.cbc.ASM.commons.Methylation.MethylStatus;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;
import static org.testng.AssertJUnit.assertTrue;

public class MappedReadTest {
	private MappedRead mappedRead;

	@BeforeMethod
	public void setUp() throws Exception {
		mappedRead = new MappedRead("test", '+', 0, "TGAACGANGA", "read1");
	}

	@Test
	public void testGetMethylStatus() throws Exception {
		// Parameterized test can be used here for large number of test cases.
		// But currently, only several simple cases added.
		assertTrue(mappedRead.getMethylStatus(0) == MethylStatus.T);
		assertTrue(mappedRead.getMethylStatus(4) == MethylStatus.C);
		assertTrue(mappedRead.getMethylStatus(7) == MethylStatus.N);
	}

	@Test
	public void testGetComplementarySequence() throws Exception {
		assertTrue("ACTTGCTNCT".equals(mappedRead.getComplementarySequence()));
	}

	@Test
	public void testToString() throws Exception {
		MappedRead read = new MappedRead("6", '+', 10, "1234567890", "test");
		assertEquals(read.toVisualizationString(5), "6\t+\t10\t20\t.....1234567890\ttest");
	}
}