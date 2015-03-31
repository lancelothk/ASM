package edu.cwru.cbc.ASM.commons;

import com.google.common.collect.Ordering;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import org.junit.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;
import static edu.cwru.cbc.ASM.commons.CommonsUtils.readReferenceGenome;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class CommonsUtilsTest {

    @Test
    public void testExtractCpGSite() throws Exception {
        List<RefCpG> refCpGList = extractCpGSite(
                "CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG", 17207805);

        assertEquals("incorrect refCpGList size", refCpGList.size(), 5);
        assertTrue("refCpGlist is unsorted!", Ordering.natural().isOrdered(refCpGList));

        // check RefCpG position
        assertEquals("incorrect RefCpG position!", refCpGList.get(0).getPos(), 17207805);
        assertEquals("incorrect RefCpG position!", refCpGList.get(1).getPos(), 17207812);
        assertEquals("incorrect RefCpG position!", refCpGList.get(2).getPos(), 17207838);
        assertEquals("incorrect RefCpG position!", refCpGList.get(3).getPos(), 17207862);
        assertEquals("incorrect RefCpG position!", refCpGList.get(4).getPos(), 17207886);
    }

    @Test
    public void testReadReferenceGenome() throws Exception {
        String testRefFileName = "testData/testRefFile.fa";
        // ref file should
        // 1. start with ">name"
        // 2. maybe single line or multiple fasta file.
        // 3. may contain upper or lower in {A,C,G,T,N}.

        // test for correct case
        RefChr refChr = readReferenceGenome(testRefFileName);
        assertEquals("chr20", refChr.getChr());
        assertEquals("ACGCAATCGNNNNNNATTGCGACGACGCGACTGNNNACGCGTAACGN", refChr.getRefString());
        assertEquals(0, refChr.getStart());
        assertEquals(46, refChr.getEnd());

        // currently ignore error cases.
    }
}