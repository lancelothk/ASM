package ASM.Execution;

import org.junit.Test;

/**
 * Created by lancelothk on 2/23/14.
 * Test checkInterval with mock data
 */
public class CheckIntervalTest {
	@Test
	public void testCheckInterval() throws Exception {
		int chrSize = 72;
		String inputFile = "ASM/testData/interval/SingleChrMappedReads";
		String outputFile = "ASM/testData/interval/SingleChrMappedReads.intervalSummary";
		CheckInterval.checkInterval(chrSize, "1", inputFile, outputFile);
	}
}
