package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.CPMR.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MappedReadLineProcessorWithSummary;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.*;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * CPMR stands for Continuous partial methylated region.
 * Main entry of the module.
 * Convert mapped reads into epigenome and split regions by continuous CpG coverage.
 */
public class CPMR {
    private static final Logger logger = Logger.getLogger(CPMR.class.getName());
    public static int MIN_CONT_COVERAGE = 4;
    public static int MIN_INTERVAL_READS = 10;
    public static int MIN_INTERVAL_CPG = 5;
    public static int outputIntervalCount = 0;

    public static void main(String[] args) throws IOException {
        String chr = "chr20";
        String ref = "hg18_" + chr + ".fa";
        String cellLine = "i90";
        String replicate = "r1";
        String homeDirectory = System.getProperty("user.home");
        String experimentPath = String.format("%s/experiments/ASM/", homeDirectory);

        String referenceGenomeFileName = experimentPath + "/data/" + ref;
        String mappedReadFileName = String.format("%s/data/%s_%s_%s", experimentPath, cellLine, replicate, chr);
        String reportPath = String.format("%s/result_%s_%s", experimentPath, cellLine, replicate);
        String outputPath = String.format("%s/result_%s_%s/intervals_%s_%d_%d_%d", experimentPath, cellLine, replicate,
                                          chr, MIN_CONT_COVERAGE, MIN_INTERVAL_CPG, MIN_INTERVAL_READS);

        setUpLogging(String.format("%s/%s-intervals_%d_%d_%d.log", reportPath, chr, MIN_CONT_COVERAGE, MIN_INTERVAL_CPG,
                                   MIN_INTERVAL_READS), Level.INFO, new SimpleFormatter());
        splitEpigenome(referenceGenomeFileName, mappedReadFileName, outputPath, reportPath, OutputFormat.MappredRead);
    }

    private static void setUpLogging(String logFileName, Level level, Formatter formatter) throws IOException {
        FileHandler outputFileHandler = new FileHandler(logFileName);
        outputFileHandler.setLevel(level);
        outputFileHandler.setFormatter(formatter);
        Logger.getLogger("").addHandler(outputFileHandler);
    }

    public static void splitEpigenome(String referenceGenomeFileName, String mappedReadFileName, String outputPath,
                                      String reportPath, OutputFormat outputFormat) throws IOException {
        long start = System.currentTimeMillis();
        File outputFile = new File(outputPath);
        if (!outputFile.exists()) {
            //noinspection ResultOfMethodCallIgnored
            outputFile.mkdirs();
        }

        RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);
        logger.info("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

        start = System.currentTimeMillis();
        refCpGList.sort((RefCpG c1, RefCpG c2) -> c1.getPos() - c2.getPos());
        logger.info("cpg sorting complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

        start = System.currentTimeMillis();
        // assign cpg id after sorted
        for (int i = 0; i < refCpGList.size(); i++) {
            refCpGList.get(i).assignIndex(i);
        }
        logger.info("cpg id assignment complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

        start = System.currentTimeMillis();
        List<MappedRead> mappedReadList = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
                                                          new MappedReadLineProcessorWithSummary(refCpGList,
                                                                                                 refChr.getRefString().length()));
        logger.info("load mappedReadList complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

        start = System.currentTimeMillis();
        List<List<RefCpG>> cpgSiteIntervalList = getIntervals(refCpGList);

        writeIntervals(outputPath, reportPath, refChr, cpgSiteIntervalList, outputFormat);
        logger.warning("Raw Interval count:\t" + cpgSiteIntervalList.size() + "");
        logger.warning("Output Interval count:\t" + outputIntervalCount + "");
        logger.info((System.currentTimeMillis() - start) / 1000.0 + "s");
    }

    private static List<List<RefCpG>> getIntervals(List<RefCpG> refCpGList) {
        // filter refCpG by coverage and partial methylation
        List<RefCpG> filteredRefCpG = refCpGList.stream().filter(
                refCpG -> refCpG.getCpGCoverage() >= MIN_CONT_COVERAGE && refCpG.hasPartialMethyl()).collect(
                Collectors.toList());

        // split regions by continuous CpG coverage
        List<List<RefCpG>> cpgSiteIntervalList = new ArrayList<>();
        boolean cont = true;
        List<RefCpG> resultRefCpGList = new ArrayList<>();
        for (int i = 0; i < filteredRefCpG.size() - 1; i++) {
            RefCpG curr = filteredRefCpG.get(i);
            if (!cont && resultRefCpGList.size() > 0) {
                cpgSiteIntervalList.add(resultRefCpGList);
                resultRefCpGList = new ArrayList<>();
            }
            resultRefCpGList.add(curr);
            RefCpG next = filteredRefCpG.get(i + 1);
            cont = curr.hasCommonRead(next);
        }
        if (cont && resultRefCpGList.size() > 0) {
            resultRefCpGList.add(resultRefCpGList.get(resultRefCpGList.size() - 1));
            cpgSiteIntervalList.add(resultRefCpGList);
        }
        return cpgSiteIntervalList;
    }

    private static void writeIntervals(String outputPath, String reportPath, RefChr refChr,
                                       List<List<RefCpG>> cpgSiteIntervalList,
                                       OutputFormat outputFormat) throws IOException {
        BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(
                String.format("%s/%s-intervalSummary_%d_%d_%d", reportPath, refChr.getChr(), MIN_CONT_COVERAGE,
                              MIN_INTERVAL_CPG, MIN_INTERVAL_READS)));
        intervalSummaryWriter.write("chr\tlength\treadCount\tCpGCount\tstartCpG\tendCpG\n");
        cpgSiteIntervalList.forEach((list) -> {
            Set<MappedRead> mappedReadSet = new HashSet<>();
            list.forEach((refCpG) -> {
                refCpG.getCpGList().forEach((CpG cpg) -> mappedReadSet.add(cpg.getMappedRead()));
                assert refCpG.getCpGList().size() >= MIN_CONT_COVERAGE; // make sure coverage is correct
            });
            // only pass high quality result for next step.
            if (mappedReadSet.size() >= MIN_INTERVAL_READS && list.size() >= MIN_INTERVAL_CPG) {
                List<MappedRead> mappedReadList = mappedReadSet.stream().sorted(
                        (m1, m2) -> m1.getId().compareTo(m2.getId())).sorted(
                        (m1, m2) -> m1.getStart() - m2.getStart()).collect(Collectors.toList());
                outputIntervalCount++;
                int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
                int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
                @SuppressWarnings("UnnecessaryLocalVariable") // since it makes clear pos and CpGPos is different.
                        int startPos = startCpGPos;
                int endPos = endCpGPos + 1;
                try {
                    intervalSummaryWriter.write(
                            String.format("%s\t%d\t%d\t%d\t%d\t%d\n", refChr.getChr(), endPos - startPos + 1,
                                          mappedReadList.size(), list.size(), startCpGPos, endCpGPos));
                    switch (outputFormat) {
                        case MappredRead:
                            Utils.writeMappedReadInInterval(outputPath, refChr, startPos, endPos, mappedReadList);
                            break;
                        case ExtEpiRead:
                            Utils.writeExtEpireadInInterval(outputPath, refChr, startPos, endPos, mappedReadList);
                            break;
                        default:
                            throw new RuntimeException("invalid output format");
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        });
        intervalSummaryWriter.close();
    }

    public enum OutputFormat {
        MappredRead, ExtEpiRead
    }
}