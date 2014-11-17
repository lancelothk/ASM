package edu.cwru.cbc.ASM.commons.DataType;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class MappedReadTest {
    private MappedRead mappedRead;

    @Before
    public void setUp() throws Exception {
        mappedRead = new MappedRead("test", '+', 0, 10, "TGAACGANGA", "read1");
    }

    @Test
    public void testGetMethylStatus() throws Exception {
        // we can use parameterized test here for large number of test cases.
        // But currently, we only use several simple cases.
        assertTrue(mappedRead.getMethylStatus(0) == MethylStatus.T);
        assertTrue(mappedRead.getMethylStatus(4) == MethylStatus.C);
        assertTrue(mappedRead.getMethylStatus(7) == MethylStatus.N);
    }

    @Test
    public void testGetComplementarySequence() throws Exception {
        assertTrue("ACTTGCTNCT".equals(mappedRead.getComplementarySequence()));
    }
}