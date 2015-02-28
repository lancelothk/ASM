package edu.cwru.cbc.ASM.commons;

import com.google.common.collect.Ordering;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import org.junit.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;
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
}