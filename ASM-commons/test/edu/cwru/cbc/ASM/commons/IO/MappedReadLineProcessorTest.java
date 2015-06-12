package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

import static org.testng.AssertJUnit.assertEquals;

public class MappedReadLineProcessorTest {
	//TODO add more unit tests

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

	@Test(dataProvider = "validCharacters")
	public void test_processLine_validCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
		List<MappedRead> mappedReadList = mlp.getResult();
		assertEquals("incorrect character!", c, mappedReadList.get(0).getSequence().charAt(0));
	}

	@Test(expectedExceptions = RuntimeException.class, dataProvider = "invalidCharacters")
	public void test_processLine_invalidCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor();
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
	}
}