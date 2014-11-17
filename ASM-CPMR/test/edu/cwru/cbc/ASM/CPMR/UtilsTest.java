package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import org.junit.Test;

import java.util.Map;

import static org.junit.Assert.assertTrue;

public class UtilsTest {

    @Test
    public void testExtractCpGSite() throws Exception {
        Map<Integer, RefCpG> refCpGMap = Utils.extractCpGSite("AACGTCGTGTCGACATCACGA");
        assertTrue(refCpGMap.containsKey(2));
        assertTrue(refCpGMap.get(2).getPos() == 2);
        assertTrue(refCpGMap.containsKey(5));
        assertTrue(refCpGMap.get(5).getPos() == 5);
        assertTrue(refCpGMap.containsKey(10));
        assertTrue(refCpGMap.get(10).getPos() == 10);
        assertTrue(refCpGMap.containsKey(18));
        assertTrue(refCpGMap.get(18).getPos() == 18);
    }
}