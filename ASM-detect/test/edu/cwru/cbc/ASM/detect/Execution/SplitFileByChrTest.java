package edu.cwru.cbc.ASM.detect.Execution;

import org.junit.Test;

/**
 * Created by ke on 2/19/14.
 */
public class SplitFileByChrTest {
	@Test
	public void testSplit() throws Exception {
		String inputFileName = "ASM/testData/split/ReadsByChr";
		String outputFilePath = "ASM/testData/split/";
		SplitFileByChr.split(inputFileName, outputFilePath);
	}
}
