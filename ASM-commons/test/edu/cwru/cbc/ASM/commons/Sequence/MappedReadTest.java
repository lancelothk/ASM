package edu.cwru.cbc.ASM.commons.Sequence;

import edu.cwru.cbc.ASM.commons.Methylation.MethylStatus;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;
import static org.testng.AssertJUnit.assertTrue;

public class MappedReadTest {
	private MappedRead plusStrandRead;
	private MappedRead minusStrandRead;

	@BeforeMethod
	public void setUp() throws Exception {
		plusStrandRead = new MappedRead("plus", '+', 0, "TGAACGANGA", "read_plus");
		minusStrandRead = new MappedRead("minus", '-', 0, "CTTTGCTCNT", "read_minus");
	}

	@Test
	public void testGetMethylStatus() throws Exception {
		assertTrue(plusStrandRead.getMethylStatus(0) == MethylStatus.T);
		assertTrue(plusStrandRead.getMethylStatus(4) == MethylStatus.C);
		assertTrue(plusStrandRead.getMethylStatus(7) == MethylStatus.N);

		assertTrue(minusStrandRead.getMethylStatus(0) == MethylStatus.T);
		assertTrue(minusStrandRead.getMethylStatus(4) == MethylStatus.C);
		assertTrue(minusStrandRead.getMethylStatus(7) == MethylStatus.N);
	}

	@Test
	public void testGetComplementarySequence() throws Exception {
		assertTrue("ACTTGCTNCT".equals(plusStrandRead.getComplementarySequence()));
		assertTrue("GAAACGAGNA".equals(minusStrandRead.getComplementarySequence()));
	}

	@Test
	public void testToString() throws Exception {
		MappedRead read = new MappedRead("6", '+', 10, "1234567890", "test");
		assertEquals(read.toVisualizationString(5), "6\t+\t10\t20\t.....1234567890\ttest");
	}
}