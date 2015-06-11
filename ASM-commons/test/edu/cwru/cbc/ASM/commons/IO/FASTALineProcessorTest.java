package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Sequence.FASTASequence;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.LinkedHashMap;

import static org.junit.Assert.assertEquals;

/**
 * Created by lancelothk on 6/10/15. Tests for FASTALineProcessor
 */
public class FASTALineProcessorTest {

	@Test(expected = RuntimeException.class)
	public void test_processLine_badStart() throws Exception {
		String badStart = "AAAAA";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(badStart);
	}

	@Test(expected = RuntimeException.class)
	public void test_processLine_emptyId() throws Exception {
		String emptyId = ">";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(emptyId);
	}

	@Test(expected = RuntimeException.class)
	public void test_processLine_invalidCharacter() throws Exception {
		String invalidCharacters = "QWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm@#$%+_=-~!^&*()/\\?><,:;\"\' ";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		for (char c : invalidCharacters.toCharArray()) {
			flp.processLine("" + c);
		}

	}

	@Test(expected = RuntimeException.class)
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

	@Test(expected = RuntimeException.class)
	public void test_getResult_duplicateID() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.processLine("ACGTN");
		flp.processLine(">test");
		flp.processLine("acgtn");
		flp.getResult();
	}

	@Test(expected = RuntimeException.class)
	public void test_processLine_missingSequence() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.processLine(">test");
	}

	@Test(expected = RuntimeException.class)
	public void test_getResult_MissingSequence() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.getResult();
	}

	@Test
	public void test_processLine_validCharacter() throws Exception {
		String validCharacter = "ACGTNacgtn.";
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test");
		flp.processLine(validCharacter);
		LinkedHashMap<String, FASTASequence> resultMap = flp.getResult();
		assertEquals("incorrect sequence!", "ACGTNacgtn.", resultMap.get("test").getSequence());
	}


	@Test
	public void test_processLine_multipleUnit() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		flp.processLine(">test1");
		flp.processLine("ACGTN");
		flp.processLine(">test2");
		flp.processLine("acgtn");
		LinkedHashMap<String, FASTASequence> resultMap = flp.getResult();
		assertEquals("incorrect number of sequence!", 2, resultMap.size());
		assertEquals("incorrect sequence!", "ACGTN", resultMap.get("test1").getSequence());
		assertEquals("incorrect sequence!", "acgtn", resultMap.get("test2").getSequence());
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
		LinkedHashMap<String, FASTASequence> resultMap = flp.getResult();
		assertEquals("incorrect number of sequence!", 2, resultMap.size());
		assertEquals("incorrect sequence!", "ACGTNACGTNACGTN", resultMap.get("test1").getSequence());
		assertEquals("incorrect sequence!", "acgtnacgtnacgtn", resultMap.get("test2").getSequence());
	}

	@Test
	public void test_getResult_emptyInput() throws Exception {
		FASTALineProcessor flp = new FASTALineProcessor();
		LinkedHashMap<String, FASTASequence> resultMap = flp.getResult();
		assertEquals("incorrect number of sequence!", 0, resultMap.size());
	}

	private <T> T getField(Object parent, Class<T> clazz,
	                       String fieldName) throws NoSuchFieldException, IllegalAccessException {
		Field f = parent.getClass().getDeclaredField(fieldName);
		f.setAccessible(true);
		return clazz.cast(f.get(parent));
	}
}