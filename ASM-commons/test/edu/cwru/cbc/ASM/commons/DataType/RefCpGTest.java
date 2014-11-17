package edu.cwru.cbc.ASM.commons.DataType;

import org.junit.Before;
import org.junit.Test;

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

        MappedRead mappedRead1 = new MappedRead("test", '+', 0, 2, "TG", "read1");
        MappedRead mappedRead2 = new MappedRead("test", '+', 0, 4, "TGTG", "read2");
        MappedRead mappedRead3 = new MappedRead("test", '+', 2, 6, "CGNG", "read3");
        MappedRead mappedRead4 = new MappedRead("test", '+', 2, 4, "CG", "read4");
        MappedRead mappedRead5 = new MappedRead("test", '+', 4, 6, "TG", "read5");

        CpG ca = new CpG(mappedRead1, rc1);
        rc1.addCpG(ca);
        ca.setMethylStatus(MethylStatus.T);
        CpG cb = new CpG(mappedRead2, rc1);
        rc1.addCpG(cb);
        cb.setMethylStatus(MethylStatus.T);
        CpG cc = new CpG(mappedRead3, rc3);
        cc.setMethylStatus(MethylStatus.N);
        rc3.addCpG(cc);
        CpG cd = new CpG(mappedRead4, rc2);
        rc2.addCpG(cd);
        cd.setMethylStatus(MethylStatus.C);
        CpG ce = new CpG(mappedRead3, rc2);
        rc2.addCpG(ce);
        ce.setMethylStatus(MethylStatus.C);
        CpG cf = new CpG(mappedRead2, rc2);
        rc2.addCpG(cf);
        cf.setMethylStatus(MethylStatus.T);
        CpG cg = new CpG(mappedRead5, rc3);
        rc3.addCpG(cg);
        cg.setMethylStatus(MethylStatus.T);
    }

    @Test
    public void testHasCommonRead() throws Exception {
        assertTrue(rc1.hasCommonRead(rc2));
        assertTrue(!rc1.hasCommonRead(rc3));
        assertTrue(rc2.hasCommonRead(rc3));
    }

    @Test
    public void testHasPartialMethyl() throws Exception {
        assertTrue(!rc1.hasPartialMethyl());
        assertTrue(rc2.hasPartialMethyl());
        assertTrue(!rc3.hasPartialMethyl());
    }
}