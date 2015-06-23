package edu.cwru.cbc.ASM.commons.sequence;

import java.util.regex.Pattern;

/**
 * Created by lancelothk on 6/20/15.
 * Class for storing IUPAC code.
 * From: http://www.bioinformatics.org/sms/iupac.html
 * Also refer to Cock, Peter JA, et al. "The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants." Nucleic acids research 38.6 (2010): 1767-1771.
 */
public class IUPACCode {
	public static final String nucleotideCode_UpperCase = "ACGTURYSWKMBDHVN.-";
	public static final String nucleotideCode_LowerCase = "acgturyswkmbdhvn.-";
	public static final String aminoAcidCode_UpperCase = "ACDEFGHIKLMNPQRSTVWY";
	public static final String aminoAcidCode_LowerCase = "acdefghiklmnpqrstvwy";
	private static final Pattern nucleotideCode_pattern = Pattern.compile(
			"[^" + Pattern.quote(nucleotideCode_UpperCase + nucleotideCode_LowerCase) + "]");
	private static final Pattern aminoAcidCode_pattern = Pattern.compile(
			"[^" + Pattern.quote(aminoAcidCode_UpperCase + aminoAcidCode_LowerCase) + "]");

	public static boolean validateNucleotideCode(String sequence) {
		return validate(nucleotideCode_pattern, sequence);
	}

	public static boolean validateAminoAcidCode(String sequence) {
		return validate(aminoAcidCode_pattern, sequence);
	}

	private static boolean validate(Pattern pattern, String sequence) {
		return !pattern.matcher(sequence).find();
	}
}
