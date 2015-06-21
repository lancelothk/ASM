package edu.cwru.cbc.ASM.commons.methylation;

import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.apache.commons.lang3.reflect.FieldUtils;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertTrue;

public class RefCpGTest {
	private RefCpG rc1;
	private RefCpG rc2;
	private RefCpG rc3;

	@BeforeMethod
	public void setUp() throws Exception {
		// ref      1       2       3
		// cpg      a(T)    d(C)    g(T)
		//          b(T)----f(T)
		//                  e(C)----c(N)

		rc1 = new RefCpG(0);
		rc2 = new RefCpG(2);
		rc3 = new RefCpG(4);

		MappedRead mappedRead1 = new MappedRead("test", '+', 0, "TG", 1);
		MappedRead mappedRead2 = new MappedRead("test", '+', 0, "TGTG", 2);
		MappedRead mappedRead3 = new MappedRead("test", '+', 2, "CGNG", 2);
		MappedRead mappedRead4 = new MappedRead("test", '+', 2, "CG", 2);
		MappedRead mappedRead5 = new MappedRead("test", '+', 4, "TG", 2);

		CpG ca = new CpG(mappedRead1, rc1, MethylStatus.T);
		rc1.addCpG(ca);
		CpG cb = new CpG(mappedRead2, rc1, MethylStatus.T);
		rc1.addCpG(cb);
		CpG cc = new CpG(mappedRead3, rc3, MethylStatus.N);
		rc3.addCpG(cc);
		CpG cd = new CpG(mappedRead4, rc2, MethylStatus.C);
		rc2.addCpG(cd);
		CpG ce = new CpG(mappedRead3, rc2, MethylStatus.C);
		rc2.addCpG(ce);
		CpG cf = new CpG(mappedRead2, rc2, MethylStatus.T);
		rc2.addCpG(cf);
		CpG cg = new CpG(mappedRead5, rc3, MethylStatus.T);
		rc3.addCpG(cg);
	}

	@Test
	public void testHasCommonRead() throws Exception {
		assertTrue(rc1.hasCommonRead(rc2, 0));
		assertTrue(!rc1.hasCommonRead(rc3, 0));
		assertTrue(rc2.hasCommonRead(rc3, 0));
	}

	@Test
	public void testHasPartialMethyl() throws Exception {
		assertTrue(!rc1.hasPartialMethyl(0.2));
		assertTrue(rc2.hasPartialMethyl(0.2));
		assertTrue(!rc3.hasPartialMethyl(0.2));
	}

	@Test
	public void testGetMajorMethylStatus() throws Exception {
		RefCpG rfCpG = new RefCpG(1);
		FieldUtils.writeField(rfCpG, "methylCount", 4, true);
		FieldUtils.writeField(rfCpG, "coveredCount", 8, true);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.E, rfCpG.getMajorMethylStatus());
		FieldUtils.writeField(rfCpG, "methylCount", 1, true);
		FieldUtils.writeField(rfCpG, "coveredCount", 8, true);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.T, rfCpG.getMajorMethylStatus());
		FieldUtils.writeField(rfCpG, "methylCount", 7, true);
		FieldUtils.writeField(rfCpG, "coveredCount", 8, true);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.C, rfCpG.getMajorMethylStatus());
	}
}