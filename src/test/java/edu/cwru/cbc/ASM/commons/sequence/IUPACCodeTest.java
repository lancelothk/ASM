package edu.cwru.cbc.ASM.commons.sequence;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashSet;

import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

public class IUPACCodeTest {

	public static final int ASCII_SIZE = 256;

	@DataProvider(name = "validNucleotideCode")
	public static Object[][] validNucleotideCode() {
		String validCharacters = IUPACCode.nucleotideCode_UpperCase + IUPACCode.nucleotideCode_LowerCase;
		Object[][] result = new Object[validCharacters.length()][];
		for (int i = 0; i < validCharacters.length(); i++) {
			result[i] = new Object[]{validCharacters.charAt(i)};
		}
		return result;
	}

	@DataProvider(name = "validAminoAcidCode")
	public static Object[][] validAminoAcidCode() {
		String validCharacters = IUPACCode.aminoAcidCode_UpperCase + IUPACCode.aminoAcidCode_LowerCase;
		Object[][] result = new Object[validCharacters.length()][];
		for (int i = 0; i < validCharacters.length(); i++) {
			result[i] = new Object[]{validCharacters.charAt(i)};
		}
		return result;
	}

	/**
	 * Use ASCII characters(256)not in validCharacter set as invalid ones.
	 */
	@DataProvider(name = "invalidNucleotideCode")
	public static Object[][] invalidNucleotideCode() {
		String validCharacters = IUPACCode.nucleotideCode_UpperCase + IUPACCode.nucleotideCode_LowerCase;
		HashSet<Character> validCharacterSet = new HashSet<>();
		for (int i = 0; i < validCharacters.length(); i++) {
			validCharacterSet.add(validCharacters.charAt(i));
		}
		Object[][] result = new Object[ASCII_SIZE - validCharacterSet.size()][];
		for (int i = 0, j = 0; i < ASCII_SIZE; i++) {
			if (!validCharacterSet.contains((char) i)) {
				result[j] = new Object[]{(char) i};
				j++;
			}
		}
		return result;
	}

	/**
	 * Use ASCII characters(256)not in validCharacter set as invalid ones.
	 */
	@DataProvider(name = "invalidAminoAcidCode")
	public static Object[][] invalidAminoAcidCode() {
		String validCharacters = IUPACCode.aminoAcidCode_UpperCase + IUPACCode.aminoAcidCode_LowerCase;
		HashSet<Character> validCharacterSet = new HashSet<>();
		for (int i = 0; i < validCharacters.length(); i++) {
			validCharacterSet.add(validCharacters.charAt(i));
		}
		Object[][] result = new Object[ASCII_SIZE - validCharacterSet.size()][];
		for (int i = 0, j = 0; i < ASCII_SIZE; i++) {
			if (!validCharacterSet.contains((char) i)) {
				result[j] = new Object[]{(char) i};
				j++;
			}
		}
		return result;
	}

	@Test(dataProvider = "validNucleotideCode")
	public void testValidateNucleotideCode(char c) throws Exception {
		assertTrue(IUPACCode.validateNucleotideCode(String.valueOf(c)));

	}

	@Test(dataProvider = "invalidNucleotideCode")
	public void testInvalidNucleotideCode(char c) throws Exception {
		assertFalse(IUPACCode.validateNucleotideCode(String.valueOf(c)));
	}

	@Test(dataProvider = "validAminoAcidCode")
	public void testValidateAminoAcidCode(char c) throws Exception {
		assertTrue(IUPACCode.validateAminoAcidCode(String.valueOf(c)));
	}

	@Test(dataProvider = "invalidAminoAcidCode")
	public void testInvalidAminoAcidCode(char c) throws Exception {
		assertFalse(IUPACCode.validateAminoAcidCode(String.valueOf(c)));
	}
}