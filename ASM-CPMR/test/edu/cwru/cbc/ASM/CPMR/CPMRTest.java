package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.ReflectionUtils.ReflectionUtils;
import org.junit.Test;

import java.io.File;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;

public class CPMRTest {

    @Test
    public void testSplitEpigenome() throws Exception {
        ReflectionUtils.setFinalStatic(CPMR.class.getField("MIN_CONT_COVERAGE"), 4);
        ReflectionUtils.setFinalStatic(CPMR.class.getField("MIN_INTERVAL_READS"), 5);
        ReflectionUtils.setFinalStatic(CPMR.class.getField("MIN_INTERVAL_CPG"), 2);

        ReflectionUtils.invokePrivateMethod(
                CPMR.class.getDeclaredMethod("setUpLogging", String.class, Level.class, Formatter.class), CPMR.class,
                "testData/test.log", Level.INFO, new Formatter() {
                    @Override
                    public String format(LogRecord record) {
                        return record.getMessage() + "\n";
                    }
                });

        CPMR.splitEpigenome("testData/i90_r1_chr20.225400_225900.test.fa", "testData/i90_r1_chr20.225500_225800.test",
                            "testData/", "testData/", CPMR.OutputFormat.MappredRead);

        junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(
                "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225518-225607.template"), new File(
                "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225518-225607"));

        junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(
                "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225722-225723.template"), new File(
                "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-225722-225723"));

        junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(
                                                         "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-intervalSummary_4_2_5.template"),
                                                 new File(
                                                         "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/chr20_225400_225900-intervalSummary_4_2_5"));

        junitx.framework.FileAssert.assertEquals("result interval file incorrect!", new File(
                "/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/test.log.template"),
                                                 new File("/home/kehu/IdeaProjects/ASM/ASM-CPMR/testData/test.log"));
    }
}