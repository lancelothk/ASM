package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Methylation.RefChr;
import org.testng.annotations.Test;

import static org.testng.AssertJUnit.assertEquals;

/**
 * Created by kehu on 6/11/15.
 * Test class for IOUtils.
 */
public class IOUtilsTest {

	@Test
	public void testReadReferenceGenome() throws Exception {
		String testRefFileName = "ASM-commons/testData/testRefFile.fa";
		// ref file should
		// 1. start with ">name"
		// 2. maybe single line or multiple fasta file.
		// 3. may contain upper or lower in {A,C,G,T,N}.

		// test for correct case
		RefChr refChr = IOUtils.readReferenceGenome(testRefFileName);
		assertEquals("incorrect reference name", "chr20", refChr.getChr());
		assertEquals("incorrecr ref string", "ACGCAATCGNNNNNNATTGCGACGACGCGACTGNNNACGCGTAACGN", refChr.getRefString());
		assertEquals("incorrect start position", 0, refChr.getStart());
		assertEquals("incorrect end position", 46, refChr.getEnd());
	}
}