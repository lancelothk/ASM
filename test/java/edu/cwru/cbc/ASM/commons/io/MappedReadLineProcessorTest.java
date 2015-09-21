package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.MappedReadFileFormat;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.testng.annotations.Test;

import java.util.List;

import static org.testng.AssertJUnit.assertEquals;

public class MappedReadLineProcessorTest {

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*invalid strand!.*")
	public void test_invalidStrand() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\tplus\t17207806\t17207874\taaaaaa\t815505";
		mlp.processLine(mappedReadStr1);
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*incompatible format detected! column number is 8*")
	public void test_invalidColumnNumber() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\teeeeee\t\t815505";
		mlp.processLine(mappedReadStr1);
	}

	@Test
	public void test_readsOrder() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\t815505";
		String mappedReadStr2 = "20\t-\t17207806\t17207874\tAAAAAA\t815506";
		String mappedReadStr3 = "20\t-\t17207806\t17207874\tAAAAAA\t815507";
		mlp.processLine(mappedReadStr1);
		mlp.processLine(mappedReadStr2);
		mlp.processLine(mappedReadStr3);
		List<MappedRead> mappedReadList = mlp.getResult();
		assertEquals("815505", mappedReadList.get(0).getId());
		assertEquals("815506", mappedReadList.get(1).getId());
		assertEquals("815507", mappedReadList.get(2).getId());
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*found duplicate mapped read!.*")
	public void test_duplicate() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\t815505";
		String mappedReadStr2 = "20\t-\t17207806\t17207874\tAAAAAA\t815505";
		mlp.processLine(mappedReadStr1);
		mlp.processLine(mappedReadStr2);
	}

	@Test
	public void test_readsFiltering() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(false, mr -> mr.getId().equals("815505"));
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\t815505";
		String mappedReadStr2 = "20\t-\t17207806\t17207874\tAAAAAA\t815506";
		mlp.processLine(mappedReadStr1);
		mlp.processLine(mappedReadStr2);
		assertEquals(1, mlp.getResult().size());
		assertEquals("815505", mlp.getResult().get(0).getId());
	}

	@Test
	public void test_pairEnd() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(true, mr -> mr.getId().startsWith("815505"));
		String mappedReadStr1 = "20\t+\t17207806\t17207827\tAAAAAA\tBBBBBB\t815505";
		mlp.processLine(mappedReadStr1);
		assertEquals("AAAAAA", mlp.getResult().get(0).getSequence());
		assertEquals(17207806, mlp.getResult().get(0).getStart());
		assertEquals("815505_1", mlp.getResult().get(0).getId());
		assertEquals("BBBBBB", mlp.getResult().get(1).getSequence());
		assertEquals(17207821, mlp.getResult().get(1).getStart());
		assertEquals("815505_2", mlp.getResult().get(1).getId());
	}

	@Test
	public void testSAMFormat() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(false, mr -> mr.getId().startsWith("815505"),
				MappedReadFileFormat.SAM);
		String mappedReadStr1 = "815505\t16\tchr20\t17207806\t255\t82M\t*\t0\t0\tAAAAAA";
		mlp.processLine(mappedReadStr1);
		assertEquals("AAAAAA", mlp.getResult().get(0).getSequence());
		assertEquals(17207806, mlp.getResult().get(0).getStart());
		assertEquals("815505", mlp.getResult().get(0).getId());
		assertEquals('-', mlp.getResult().get(0).getStrand());
	}
}