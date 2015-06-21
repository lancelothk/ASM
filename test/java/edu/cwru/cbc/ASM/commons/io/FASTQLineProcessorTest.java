package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.sequence.FASTQSequence;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.LinkedHashSet;

import static org.testng.Assert.assertEquals;

/**
 * Created by kehu on 6/11/15. Tests for FASTALineProcessor
 */
public class FASTQLineProcessorTest {

	@DataProvider(name = "invalidQualityScores")
	public static Object[][] invalidQualityScores() {
		Object[][] result = new Object[40][];
		for (int i = 0; i < 33; i++) {
			result[i] = new Object[]{i};
		}
		for (int i = 33; i < 40; i++) {
			result[i] = new Object[]{i + 94};
		}
		return result;
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*unexpected line.*")
	public void test_processLine_missingIDInTheBeginnning() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("ACGT");
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*unexpected line.*")
	public void test_processLine_missingID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("+test");
		flp.processLine("ffff");
		flp.processLine("ACGT");
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*duplicate id detected.*")
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


	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*invalid character in sequence.*")
	public void test_processLine_missingSequence() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("+test");
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*unexpected line.*")
	public void test_processLine_missingPlus() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("ffff");
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*length of quality String does not match length of sequence.*")
	public void test_processLine_missingQuality() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.processLine("+test");
		flp.processLine("@test");
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*ID in plus line is not same to ID in ID line.*")
	public void test_processLine_invalidPlusID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("acgt");
		flp.processLine("+something");
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*invalid quality score character.*",
			dataProvider = "invalidQualityScores")
	public void test_processLine_invalidQualityScore(final int c) throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("a");
		flp.processLine("+test");
		flp.processLine("" + (char) c);
	}


	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*empty ID line.*")
	public void test_processLine_emptyID() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@");
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*missing sequence line in the end of file!.*")
	public void test_getResult_missingSequence() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.getResult();
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*missing plus line in the end of file!.*")
	public void test_getResult_missingPlus() throws Exception {
		FASTQLineProcessor flp = new FASTQLineProcessor();
		flp.processLine("@test");
		flp.processLine("ACGT");
		flp.getResult();
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*missing quality line in the end of file!.*")
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