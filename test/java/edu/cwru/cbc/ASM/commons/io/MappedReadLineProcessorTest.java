package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

import static org.testng.AssertJUnit.assertEquals;

public class MappedReadLineProcessorTest {
	@DataProvider(name = "invalidCharacters")
	public static Object[][] invalidCharacters() {
		String invalidCharacters = "acgtnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm@#$%+_=-~!^&*()/\\?><,:;\"\' ";
		Object[][] result = new Object[invalidCharacters.length()][];
		for (int i = 0; i < invalidCharacters.length(); i++) {
			result[i] = new Object[]{invalidCharacters.charAt(i)};
		}
		return result;
	}

	@DataProvider(name = "validCharacters")
	public static Object[][] validCharacters() {
		String validCharacters = "ACGTN.";
		Object[][] result = new Object[validCharacters.length()][];
		for (int i = 0; i < validCharacters.length(); i++) {
			result[i] = new Object[]{validCharacters.charAt(i)};
		}
		return result;
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*invalid strand!.*")
	public void test_invalidStrand() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\tplus\t17207806\t17207874\taaaaaa\t815505";
		mlp.processLine(mappedReadStr1);
	}

	@Test(dataProvider = "validCharacters")
	public void test_processLine_validCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
		List<MappedRead> mappedReadList = mlp.getResult();
		assertEquals("incorrect character!", c, mappedReadList.get(0).getSequence().charAt(0));
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*invalid character in sequence!.*", dataProvider = "invalidCharacters")
	public void test_processLine_invalidCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*columns is not correct for mapped read format.*")
	public void test_invalidColumnNumber() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\teeeeee\t\t815505";
		mlp.processLine(mappedReadStr1);
	}

	@Test
	public void test_sevenColumnNumber() throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\teeeeee\t815505";
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
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(mr -> mr.getId().equals("815505"));
		String mappedReadStr1 = "20\t-\t17207806\t17207874\tAAAAAA\t815505";
		String mappedReadStr2 = "20\t-\t17207806\t17207874\tAAAAAA\t815506";
		mlp.processLine(mappedReadStr1);
		mlp.processLine(mappedReadStr2);
		assertEquals(1, mlp.getResult().size());
		assertEquals("815505", mlp.getResult().get(0).getId());
	}
}