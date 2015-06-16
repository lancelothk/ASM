package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Sequence.FASTQSequence;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.LinkedHashSet;

import static org.testng.Assert.assertEquals;

/**
 * Created by kehu on 6/11/15. Tests for FASTALineProcessor
 */
public class FASTQLineProcessorTest {

	@DataProvider(name = "invalidCharacters")
	public static Object[][] invalidCharacters() {
		String invalidCharacters = "QWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm@#$%+_=-~!^&*()/\\?><,:;\"\' ";
		Object[][] result = new Object[invalidCharacters.length()][];
		for (int i = 0; i < invalidCharacters.length(); i++) {
			result[i] = new Object[]{invalidCharacters.charAt(i)};
		}
		return result;
	}

	@DataProvider(name = "invalidQualityScores")
	public static Object[][] invalidQualityScores() {
		Object[][] result = new Object[33][];
		for (int i = 0; i < 33; i++) {
			result[i] = new Object[]{i};
		}
		return result;
	}

	@DataProvider(name = "validCharacters")
	public static Object[][] validCharacters() {
		String validCharacters = "acgtnACGTN.";
		Object[][] result = new Object[validCharacters.length()][];
		for (int i = 0; i < validCharacters.length(); i++) {
			result[i] = new Object[]{validCharacters.charAt(i)};
		}
		return result;
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_missingIDInTheBeginnning() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("ACGT");
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_missingID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("+test");
		flp.processLine("ffff");
		flp.processLine("ACGT");
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_DuplicateID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGTC");
		flp.processLine("+test");
		flp.processLine("ffffe");
		flp.processLine("@test");
		flp.processLine("ACGTA");
		flp.processLine("+test");
		flp.processLine("ffffe");
	}


	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_missingSequence() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("+test");
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_missingPlus() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("ffff");
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_missingQuality() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("+test");
		flp.processLine("@test2");
	}

	@Test(expectedExceptions = RuntimeException.class, dataProvider = "invalidCharacters")
	public void test_processLine_invalidSequenceCharacter(final char c) throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@" + c);
		flp.processLine("" + c);
		flp.processLine("+" + c);
		flp.processLine("f");
	}


	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_invalidPlusID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("acgt");
		flp.processLine("+something");
	}

	@Test(expectedExceptions = RuntimeException.class, dataProvider = "invalidQualityScores")
	public void test_processLine_invalidQualityScore(final int c) throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("a");
		flp.processLine("+test");
		flp.processLine("" + (char) c);
	}


	@Test(expectedExceptions = RuntimeException.class)
	public void test_processLine_emptyID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@");
	}

	@Test(dataProvider = "validCharacters")
	public void test_processLine_validCharacter(char c) throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@" + c);
		flp.processLine("" + c);
		flp.processLine("+" + c);
		flp.processLine("f");
		LinkedHashSet<FASTQSequence> resultSet = flp.getResult();
		assertEquals(resultSet.iterator().next().getSequence().charAt(0), c, "incorrect character!");
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_getResult_missingSequence() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.getResult();
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_getResult_missingPlus() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.getResult();
	}

	@Test(expectedExceptions = RuntimeException.class)
	public void test_getResult_missingQuality() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("+test");
		flp.getResult();
	}

	@Test
	public void test_getResult_emptyInput() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		LinkedHashSet<FASTQSequence> resultSet = flp.getResult();
		assertEquals(resultSet.size(), 0, "incorrect number of sequence!");
	}
}