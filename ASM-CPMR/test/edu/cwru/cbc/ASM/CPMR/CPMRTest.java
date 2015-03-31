package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import org.junit.Test;

public class CPMRTest {

	@Test
	public void testSplitEpigenome() throws Exception {
		String homeDirectory = System.getProperty("user.home");
		ReflectionUtils.setFinalStaticField(CPMR.class.getField("MIN_CONT_COVERAGE"), 4);
		ReflectionUtils.setFinalStaticField(CPMR.class.getField("MIN_INTERVAL_READS"), 5);
		ReflectionUtils.setFinalStaticField(CPMR.class.getField("MIN_INTERVAL_CPG"), 2);

		//		ReflectionUtils.invokePrivateMethod(
		//				CPMR.class.getDeclaredMethod("setUpLogging", String.class, Level.class, Formatter.class), CPMR.class,
		//				homeDirectory + "/IdeaProjects/ASM/ASM-CPMR/testData/test.log", Level.WARNING, new Formatter() {
		//					@Override
		//					public String format(LogRecord record) {
		//						return record.getMessage() + "\n";
		//					}
		//				});

		//TODO fix this unit test

		//		CPMR.splitEpigenome(homeDirectory + "/IdeaProjects/ASM/ASM-CPMR/testData/i90_r1_chr20.225400_225900.test.fa",
		//							homeDirectory + "/IdeaProjects/ASM/ASM-CPMR/testData/i90_r1_chr20.225500_225800.test",
		//							homeDirectory + "/IdeaProjects/ASM/ASM-CPMR/testData/",
		//							homeDirectory + "/IdeaProjects/ASM/ASM-CPMR/testData/", CPMR.OutputFormat.);

		//		junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(homeDirectory +
		//																									 "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225518-225607.template"),
		//												 new File(homeDirectory +
		//																  "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225518-225607"));
		//
		//		junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(homeDirectory +
		//																									 "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225722-225723.template"),
		//												 new File(homeDirectory +
		//																  "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225722-225723"));
		//
		//		junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(homeDirectory +
		//																									 "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-intervalSummary_4_2_5.template"),
		//												 new File(homeDirectory +
		//																  "/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-intervalSummary_4_2_5"));
		//
		//		junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(homeDirectory +
		//																									 "/IdeaProjects/ASM/ASM-CPMR/testData/test.log.template"),
		//												 new File(homeDirectory +
		//																  "/IdeaProjects/ASM/ASM-CPMR/testData/test.log"));
	}
}