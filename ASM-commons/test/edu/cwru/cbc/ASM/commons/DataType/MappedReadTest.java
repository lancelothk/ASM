package edu.cwru.cbc.ASM.commons.DataType;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class MappedReadTest {
    private MappedRead mappedRead;

    @Before
    public void setUp() throws Exception {
        mappedRead = new MappedRead("test", '+', 0, 10, "TGAACGANGA", "read1");
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
        MappedRead read = new MappedRead("6", '+', 10, 20, "1234567890", "test");
        assertEquals(read.toVisualizationString(5), "6\t+\t10\t20\t.....1234567890\ttest");
        System.out.println(read.toVisualizationString(5));
    }
}