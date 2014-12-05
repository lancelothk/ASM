package edu.cwru.cbc.ASM.CPMR;

import org.junit.Test;

public class CPMRTest {

    @Test
    public void testSplitEpigenome() throws Exception {
        CPMR.splitEpigenome("testData/i90_r1_chr20.225400_225900.test.fa", "testData/i90_r1_chr20.225500_225800.test",
                            "testData/", "testData/", CPMR.OutputFormat.MappredRead);
    }
}