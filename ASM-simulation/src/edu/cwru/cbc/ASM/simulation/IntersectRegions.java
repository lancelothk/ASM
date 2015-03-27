package edu.cwru.cbc.ASM.simulation;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.DataType.GenomicInterval;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegions {
    public static void main(String[] args) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("t", true, "target file");
        options.addOption("r", true, "actual result file");
        options.addOption("e", true, "expected result file");
        options.addOption("o", true, "output file");

        CommandLineParser parser = new BasicParser();
        CommandLine cmd = parser.parse(options, args);

        String targetFileName = cmd.getOptionValue("t");
        String actualResultFileName = cmd.getOptionValue("r");
        String expectedResultFileName = cmd.getOptionValue("e");
        String outputFileName = cmd.getOptionValue("o");

        execution(targetFileName, actualResultFileName, expectedResultFileName, outputFileName);
    }

    public static void execution(String targetFileName, String actualResultFileName, String expectedResultFileName,
                                 String outputFileName) throws IOException {
        System.out.println("actualResultFile:\t" + actualResultFileName);
        twoWayIntersection(targetFileName, actualResultFileName, expectedResultFileName);
        fourWayIntersection(actualResultFileName, expectedResultFileName, outputFileName);
        System.out.println();
    }


    public static void calcOverlappedLenth(String targetFileName, String actualResultFileName,
                                           String label) throws IOException {
        int overlappedLength = 0, tpLength = 0;
        List<GenomicInterval> targetRegions = CommonsUtils.readBedRegions(targetFileName);
        List<GenomicInterval> resultRegions = CommonsUtils.readBedRegions(actualResultFileName + label);

        for (GenomicInterval resultRegion : resultRegions) {
            tpLength += resultRegion.length();
        }

        for (GenomicInterval targetRegion : targetRegions) {
            for (GenomicInterval resultRegion : resultRegions) {
                if (targetRegion.getStart() <= resultRegion.getEnd() &&
                        targetRegion.getEnd() >= resultRegion.getStart()) {
                    targetRegion.addIntersectedRegion(resultRegion);
                    resultRegion.addIntersectedRegion(targetRegion);
                    int right = Math.min(targetRegion.getEnd(), resultRegion.getEnd());
                    int left = Math.max(targetRegion.getStart(), resultRegion.getStart());
                    int length = right - left + 1;
                    System.out.println(resultRegion.getStart() + "\t" + resultRegion.getEnd() +
                                               "\toverlapped length/" + label + " length:\t" + length + "\t" +
                                               resultRegion.length());
                    overlappedLength += length;
                }
            }
        }
        System.out.println("overlapped length:\t" + overlappedLength);
        System.out.println(label + " length:\t" + tpLength);
    }

    private static void twoWayIntersection(String targetFileName, String actualResultFileName,
                                           String expectedResultFileName) throws IOException {
        List<GenomicInterval> targetRegions = CommonsUtils.readBedRegions(targetFileName);
        List<GenomicInterval> resultRegions = CommonsUtils.readBedRegions(actualResultFileName);

        BufferedWriter writer = new BufferedWriter(new FileWriter(expectedResultFileName));

        for (GenomicInterval targetRegion : targetRegions) {
            for (GenomicInterval resultRegion : resultRegions) {
                if (targetRegion.getStart() <= resultRegion.getEnd() &&
                        targetRegion.getEnd() >= resultRegion.getStart()) {
                    targetRegion.addIntersectedRegion(resultRegion);
                    resultRegion.addIntersectedRegion(targetRegion);
                }
            }
        }
        for (GenomicInterval resultRegion : resultRegions) {
            if (resultRegion.isIntersected()) {
                writer.write(String.format("%s\t", resultRegion.toBedString()));
                for (GenomicInterval intersectedRegion : resultRegion.getIntersectedRegions()) {
                    writer.write(String.format("%s\t", intersectedRegion.toBedString()));
                }
                writer.write("+\n");
            }
        }

        for (GenomicInterval targetRegion : targetRegions) {
            if (!targetRegion.isIntersected()) {
                writer.write(targetRegion.toBedString() + "\t*\n");
            }
        }

        for (GenomicInterval resultRegion : resultRegions) {
            if (!resultRegion.isIntersected()) {
                writer.write(resultRegion.toBedString() + "\t-\n");
            }
        }

        writer.close();
    }

    private static void fourWayIntersection(String actualResultFileName, String expectedResultFileName,
                                            String outputFileName) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
        int tp = 0, fp = 0, fn = 0, tn = 0;
        List<GenomicInterval> expectedRegions = CommonsUtils.readBedRegions(expectedResultFileName, true);
        List<GenomicInterval> actualRegions = CommonsUtils.readBedRegions(actualResultFileName, true);
        for (GenomicInterval expectedRegion : expectedRegions) {
            for (GenomicInterval actualRegion : actualRegions) {
                if (expectedRegion.getStart() == actualRegion.getStart() &&
                        expectedRegion.getEnd() == actualRegion.getEnd()) {
                    if (expectedRegion.isPositive() && actualRegion.isPositive()) {
                        tp++;
                        writer.write(expectedRegion.toBedString() + "\ttp\n");
                    }
                    if (expectedRegion.isPositive() && !actualRegion.isPositive()) {
                        fn++;
                        writer.write(expectedRegion.toBedString() + "\tfn\n");
                    }
                    if (!expectedRegion.isPositive() && actualRegion.isPositive()) {
                        fp++;
                        writer.write(expectedRegion.toBedString() + "\tfp\n");
                    }
                    if (!expectedRegion.isPositive() && !actualRegion.isPositive()) {
                        tn++;
                        writer.write(expectedRegion.toBedString() + "\ttn\n");
                    }
                }
            }
        }
        writer.close();
        System.out.printf("truePositive:%d\ttrueNegative:%d\n", tp, tn);
        System.out.printf("falseNegative:%d\tfalsePositive:%d\n", fn, fp);
        System.out.printf("accuracy:%.4f\n", (tp + tn) / (double) (tp + tn + fn + fp));
        System.out.printf("Precision:%.4f\n", tp / (double) (tp + fp));
        System.out.printf("Sensitivity:%.4f\n", tp / (double) (tp + fn));
        System.out.printf("Fall-out:%.4f\n", fp / (double) (fp + tn));
        System.out.printf("FDR:%.4f\n", fp / (double) (fp + tp));
        System.out.printf("F1 score:%.4f\n", 2 * tp / (double) (2 * tp + fp + fn));
        System.out.printf("MCC:%.4f\n", (tp * tn - fp * fn) / (double) ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)));
    }
}
