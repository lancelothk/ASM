package edu.cwru.cbc.ASM.commons;

import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import org.junit.Test;

import java.util.List;

public class UtilsTest {

    @Test
    public void testExtractCpGSite() throws Exception {
        List<RefCpG> refCpGList = Utils.extractCpGSite("AACGTCGTGTCGACATCACGA", 0);
        //TODO
//        assertTrue(refCpGList.containsKey(2));
//        assertTrue(refCpGList.get(2).getPos() == 2);
//        assertTrue(refCpGList.containsKey(5));
//        assertTrue(refCpGList.get(5).getPos() == 5);
//        assertTrue(refCpGList.containsKey(10));
//        assertTrue(refCpGList.get(10).getPos() == 10);
//        assertTrue(refCpGList.containsKey(18));
//        assertTrue(refCpGList.get(18).getPos() == 18);
    }
}