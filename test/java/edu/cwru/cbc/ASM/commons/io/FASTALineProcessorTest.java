package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.sequence.FASTASequence;
import org.testng.annotations.Test;

import java.util.LinkedHashSet;

import static org.testng.Assert.assertEquals;


/**
 * Created by lancelothk on 6/10/15. Tests for FASTALineProcessor
 */
public class FASTALineProcessorTest {

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*missing '>' in the beginning of file.*")
	public void test_processLine_missingIdInBeginning() throws Exception {
		String missingIdInBeginning = "AAAAA";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(missingIdInBeginning);
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*empty id!.*")
	public void test_processLine_emptyId() throws Exception {
		String emptyId = ">";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(emptyId);
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*duplicate id.*")
	public void test_processLine_duplicateID() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		// Since the last sequence unit is put into map only when getResult called, processLine cannot detect duplicate
		// in last seuqnce unit. So use 3 units here.
		flp.processLine(">test");
		flp.processLine("ACGTN");
		flp.processLine(">test");
		flp.processLine("acgtn");
		flp.processLine(">test");
		flp.processLine("ACGTN");
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*missing sequence line.*")
	public void test_processLine_missingSequence() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.processLine(">test");
	}

	@Test(expectedExceptions = RuntimeException.class, expectedExceptionsMessageRegExp = ".*missing sequence line.*")
	public void test_getResult_MissingSequence() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.getResult();
	}

	@Test
	public void test_processLine_multipleUnit() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test1");
		flp.processLine("ACGTN");
		flp.processLine(">test2");
		flp.processLine("acgtn");
		LinkedHashSet<FASTASequence> resultSet = flp.getResult();
		FASTASequence[] result = resultSet.toArray(new FASTASequence[2]);
		assertEquals(2, resultSet.size(), "incorrect number of sequence!");
		assertEquals("ACGTN", result[0].getSequence(), "incorrect sequence!");
		assertEquals("acgtn", result[1].getSequence(), "incorrect sequence!");
	}

	@Test
	public void test_processLine_multipleLineSequence() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test1");
		flp.processLine("ACGTN");
		flp.processLine("ACGTN");
		flp.processLine("ACGTN");
		flp.processLine(">test2");
		flp.processLine("acgtn");
		flp.processLine("acgtn");
		flp.processLine("acgtn");
		LinkedHashSet<FASTASequence> resultSet = flp.getResult();
		assertEquals(2, resultSet.size(), "incorrect number of sequence!");
		FASTASequence[] result = resultSet.toArray(new FASTASequence[2]);
		assertEquals("ACGTNACGTNACGTN", result[0].getSequence(), "incorrect sequence!");
		assertEquals("acgtnacgtnacgtn", result[1].getSequence(), "incorrect sequence!");
	}

	@Test
	public void test_getResult_emptyInput() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		LinkedHashSet<FASTASequence> resultSet = flp.getResult();
		assertEquals(0, resultSet.size(), "incorrect number of sequence!");
	}

}