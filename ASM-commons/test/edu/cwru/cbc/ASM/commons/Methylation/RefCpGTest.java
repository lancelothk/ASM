package edu.cwru.cbc.ASM.commons.Methylation;

import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class RefCpGTest {
	private RefCpG rc1;
	private RefCpG rc2;
	private RefCpG rc3;

	@Before
	public void setUp() throws Exception {
		// ref      1       2       3
		// cpg      a(T)    d(C)    g(T)
		//          b(T)----f(T)
		//                  e(C)----c(N)

		rc1 = new RefCpG(0);
		rc2 = new RefCpG(2);
		rc3 = new RefCpG(4);

		MappedRead mappedRead1 = new MappedRead("test", '+', 0, "TG", "read1");
		MappedRead mappedRead2 = new MappedRead("test", '+', 0, "TGTG", "read2");
		MappedRead mappedRead3 = new MappedRead("test", '+', 2, "CGNG", "read3");
		MappedRead mappedRead4 = new MappedRead("test", '+', 2, "CG", "read4");
		MappedRead mappedRead5 = new MappedRead("test", '+', 4, "TG", "read5");

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
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), rfCpG, 4);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), rfCpG, 8);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.E, rfCpG.getMajorMethylStatus());
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), rfCpG, 1);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), rfCpG, 8);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.T, rfCpG.getMajorMethylStatus());
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), rfCpG, 7);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), rfCpG, 8);
		assertEquals("incorrect MajorMethylStatus!", MethylStatus.C, rfCpG.getMajorMethylStatus());
	}
}