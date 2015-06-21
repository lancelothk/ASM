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

	public static boolean validateNucleotideCode(String sequence) {
		return validate(sequence, nucleotideCode_UpperCase + nucleotideCode_LowerCase);
	}

	public static boolean validateAminoAcidCode(String sequence) {
		return validate(sequence, aminoAcidCode_UpperCase + aminoAcidCode_LowerCase);
	}

	private static boolean validate(String sequence, String code) {
		return !Pattern.compile("[^" + Pattern.quote(code) + "]").matcher(sequence).find();
	}
}
