package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.genomicInterval.ImmutableGenomicInterval;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.testng.AssertJUnit.assertEquals;

public class CPMRTest {
	RefChr refChr;

	@BeforeMethod
	public void setUp() throws Exception {
		refChr = mock(RefChr.class);
		when(refChr.getChr()).thenReturn("chr20");
		when(refChr.getRefString()).thenReturn("CGCGCGCGCG");
	}

	@Test
	public void testGetGenomicIntervals_oneRegion() throws Exception {
		List<RefCpG> refCpGList = new ArrayList<>();
		refCpGList.add(new RefCpG(0));
		refCpGList.add(new RefCpG(2));
		refCpGList.add(new RefCpG(4));
		refCpGList.add(new RefCpG(6));
		refCpGList.add(new RefCpG(8));
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		mlp.processLine("20\t+\t0\t3\tCGCG\t0");
		mlp.processLine("20\t+\t2\t5\tCGCG\t1");
		mlp.processLine("20\t+\t4\t7\tTGCG\t2");
		mlp.processLine("20\t+\t2\t5\tTGCG\t3");
		mlp.processLine("20\t+\t4\t7\tTGTG\t4");
		mlp.processLine("20\t+\t6\t9\tTGTG\t5");
		List<MappedRead> mappedReadList = mlp.getResult();
		Map<Integer, RefCpG> refCpGMap = refCpGList.stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refCpGMap));

		CPMR cpmr = new CPMR(refCpGList, refChr, 2, 3, 4, 0.2);
		List<ImmutableGenomicInterval> resultIntervalList = cpmr.getGenomicIntervals();
		assertEquals(1, cpmr.getRawIntervalCount());
		assertEquals(1, resultIntervalList.size());
		assertEquals(6, resultIntervalList.get(0).getMappedReadList().size());
		assertEquals(3, resultIntervalList.get(0).getRefCpGList().size());
	}

	@Test
	public void testGetGenomicIntervals_twoRegions() throws Exception {
		List<RefCpG> refCpGList = new ArrayList<>();
		refCpGList.add(new RefCpG(0));
		refCpGList.add(new RefCpG(2));
		refCpGList.add(new RefCpG(4));
		refCpGList.add(new RefCpG(6));
		refCpGList.add(new RefCpG(8));
		refCpGList.add(new RefCpG(10));
		refCpGList.add(new RefCpG(12));
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		mlp.processLine("20\t+\t2\t5\tCGCG\t1");
		mlp.processLine("20\t+\t2\t5\tTGTG\t2");
		mlp.processLine("20\t+\t8\t11\tTGTG\t3");
		mlp.processLine("20\t+\t8\t11\tCGCG\t4");
		mlp.processLine("20\t+\t0\t3\tCGCG\t5");
		List<MappedRead> mappedReadList = mlp.getResult();
		Map<Integer, RefCpG> refCpGMap = refCpGList.stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refCpGMap));

		CPMR cpmr = new CPMR(refCpGList, refChr, 2, 2, 2, 0.2);
		List<ImmutableGenomicInterval> resultIntervalList = cpmr.getGenomicIntervals();
		assertEquals(2, cpmr.getRawIntervalCount());
		assertEquals(2, resultIntervalList.size());
		assertEquals(3, resultIntervalList.get(0).getMappedReadList().size());
		assertEquals(2, resultIntervalList.get(0).getRefCpGList().size());
		assertEquals(2, resultIntervalList.get(1).getMappedReadList().size());
		assertEquals(2, resultIntervalList.get(1).getRefCpGList().size());
	}
}