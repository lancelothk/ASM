package edu.cwru.cbc.ASM.commons;

import com.google.common.collect.Ordering;
import edu.cwru.cbc.ASM.commons.DataType.GenomicInterval;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.*;
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
		assertEquals("incorrect reference name", "chr20", refChr.getChr());
		assertEquals("incorrecr ref string", "ACGCAATCGNNNNNNATTGCGACGACGCGACTGNNNACGCGTAACGN", refChr.getRefString());
		assertEquals("incorrect start position", 0, refChr.getStart());
		assertEquals("incorrect end position", 46, refChr.getEnd());

		// TODO add error cases.
	}

	@Test
	public void testReadBedRegions() throws Exception {
		String bedWithLabel = "testData/bedWithLabel.bed";

		Map<String, List<GenomicInterval>> bedRegionsWithLabel = readBedRegions(bedWithLabel, true);

		assertEquals("bed region size incorrect!", 11, bedRegionsWithLabel.get("chr20").size());
		assertEquals("incorrect label!", true, bedRegionsWithLabel.get("chr20").get(0).isPositive());


		// TODO add error cases.
	}
}