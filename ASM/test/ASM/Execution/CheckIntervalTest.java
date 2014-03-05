package ASM.Execution;

import ASM.DataType.ChrCoverageSummary;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Created by lancelothk on 2/23/14.
 * Test checkInterval with mock data
 */
public class CheckIntervalTest {
	@Test
	public void testCheckInterval() throws Exception {
		int chrSize = 72;
		String referenceFile = "ASM/testData/interval/testRef";
		String inputFile = "ASM/testData/interval/SingleChrMappedReads";
		String outputFile = "ASM/testData/interval/SingleChrMappedReads.intervalSummary";
		CheckInterval.checkInterval(chrSize, "1", referenceFile, inputFile, outputFile);
	}

	@Test
	public void testCpGCount() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException, IOException {
		Method countCpG = ChrCoverageSummary.class.getDeclaredMethod("countCpG", String.class, Integer.TYPE, Integer.TYPE);
		countCpG.setAccessible(true);
		ChrCoverageSummary emptyChrCoverageSummary = new ChrCoverageSummary(0);
		String reference = "acGgCtCggccCc";
		Assert.assertEquals(1, countCpG.invoke(emptyChrCoverageSummary, reference, 1, 5));
		Assert.assertEquals(0, countCpG.invoke(emptyChrCoverageSummary, reference, 8, 12));
		Assert.assertEquals(0, countCpG.invoke(emptyChrCoverageSummary, reference, 8, 13));
		Assert.assertEquals(0, countCpG.invoke(emptyChrCoverageSummary, reference, 8, 14));
		Assert.assertEquals(2, countCpG.invoke(emptyChrCoverageSummary, reference, 1, 12));
	}
}
