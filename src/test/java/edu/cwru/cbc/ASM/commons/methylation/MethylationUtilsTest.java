package edu.cwru.cbc.ASM.commons.methylation;

import com.google.common.collect.Ordering;
import org.testng.annotations.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertTrue;

public class MethylationUtilsTest {

	@Test
	public void testExtractCpGSite() throws Exception {
		// covered cases: 1. CG in the beginning; 2.CG in the end; 3. CG in the middle.
		List<RefCpG> refCpGList = extractCpGSite(
				"CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG", 17207804);

		assertEquals("incorrect refCpGList size", refCpGList.size(), 5);
		assertTrue("refCpGlist is unsorted!", Ordering.natural().isOrdered(refCpGList));

		// check RefCpG position
		assertEquals("incorrect RefCpG position!", refCpGList.get(0).getPos(), 17207804);
		assertEquals("incorrect RefCpG position!", refCpGList.get(1).getPos(), 17207811);
		assertEquals("incorrect RefCpG position!", refCpGList.get(2).getPos(), 17207837);
		assertEquals("incorrect RefCpG position!", refCpGList.get(3).getPos(), 17207861);
		assertEquals("incorrect RefCpG position!", refCpGList.get(4).getPos(), 17207885);
	}

	@Test
	public void test_ExtractCpGSite_shortRef() throws Exception {
		List<RefCpG> refCpGList = extractCpGSite("CG", 17207804);
		assertEquals("incorrect refCpGList size", refCpGList.size(), 1);
		assertEquals("incorrect RefCpG position!", refCpGList.get(0).getPos(), 17207804);

		refCpGList = extractCpGSite("G", 17207805);
		assertEquals("incorrect refCpGList size", refCpGList.size(), 0);
	}
}