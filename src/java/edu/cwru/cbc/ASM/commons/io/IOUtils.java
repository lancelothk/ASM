package edu.cwru.cbc.ASM.commons.io;

import com.google.common.io.CharStreams;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.sequence.FASTASequence;
import org.apache.commons.lang3.StringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashSet;

/**
 * Created by kehu on 6/11/15.
 * Utils class for IO related methods.
 */
public class IOUtils {

	/**
	 * Read Reference genome into RefChr type.
	 * The ref string is concerted to upper case.
	 */
	public static RefChr readReferenceGenome(String inputFileName) throws IOException {
		File refFile = new File(inputFileName);
		int initCapacity;
		if (refFile.length() <= Integer.MAX_VALUE) {
			initCapacity = (int) refFile.length();
		} else {
			initCapacity = Integer.MAX_VALUE;
		}
		LinkedHashSet<FASTASequence> sequences = CharStreams.readLines(new FileReader(refFile), new FASTALineProcessor(initCapacity));
		if (sequences.size() == 0) {
			throw new RuntimeException("No sequence read from reference file!");
		} else if (sequences.size() > 1) {
			throw new RuntimeException("Reference file contains multiple sequences(ids)!");
		}
		FASTASequence fastaSequence = sequences.iterator().next();
		return new RefChr(fastaSequence.getId(), fastaSequence.getSequence().toUpperCase());
	}

	public static String readRefFromIntervalReadsFile(File inputFile) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFile));
		String line = bufferedReader.readLine();
		String[] items = line.split("\t");
		if (items.length != 2 || !items[0].equals("ref")) {
			throw new RuntimeException(
					"invalid reference line in interval read file!\t" + line + "\t" + inputFile.getName());
		} else {
			// item[1](reference string) should be in uppercase and without space
			assert StringUtils.isAllUpperCase(items[1]);
			assert !items[1].contains(" ");
			return items[1];
		}
	}
}
