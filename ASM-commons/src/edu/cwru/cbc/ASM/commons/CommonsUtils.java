package edu.cwru.cbc.ASM.commons;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by kehu on 11/13/14.
 * CommonsUtils for methods shared by modules under ASM project.
 */
public class CommonsUtils {

	/**
	 * Extract reference CpG from reference string.
	 *
	 * @param reference should be upperCase string and without space in it.
	 * @return returned List is sorted by position.
	 */
	public static List<RefCpG> extractCpGSite(String reference, int initPos) {
		// initialize refCpGList with 1/20 of reference size, since the probability of CG occurrence should be 1/16.
		// Actual probability is higher than 1/16 from current observation. In CPG dense region, it should be much higher.
		List<RefCpG> reFCpGList = new ArrayList<>(reference.length() / 20);
		for (int i = 0; i < reference.length() - 1; i++) {
			if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
				reFCpGList.add(new RefCpG(i + initPos));
			}
		}
		return reFCpGList;
	}

	/**
	 * Read Reference genome into RefChr type.
	 * The ref string is concerted to upper case.
	 *
	 * @throws IOException
	 */
	public static RefChr readReferenceGenome(String inputFileName) throws IOException {
		return Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor<RefChr>() {
			private String chr;
			private StringBuilder referenceBuilder = new StringBuilder();

			@Override
			public boolean processLine(String s) throws IOException {
				if (s.startsWith(">")) {
					chr = s.substring(1, s.length());
				} else {
					if (Pattern.compile("[^ACGTNacgtn]").matcher(s).find()) {
						throw new RuntimeException("invalid character in sequence! only acgtn/ACGTN are allowed!");
					}
					referenceBuilder.append(s.replaceAll(" ", "").toUpperCase());
				}
				return true;
			}

			@Override
			public RefChr getResult() {
				if (chr == null) {
					throw new RuntimeException("didn't find reference name!");
				}
				return new RefChr(chr, referenceBuilder.toString());
			}
		});
	}
}
