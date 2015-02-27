package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.DataType.GenomicRegion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.Utils.readBedRegions;

/**
 * Created by kehu on 2/16/15.
 * check intersection of two region groups
 */
public class IntersectRegions {
    public static void main(String[] args) throws IOException {
        String currUserHome = System.getProperty("user.home");
        execution(currUserHome, 1, 0);
        execution(currUserHome, 0.8, 0.2);
        execution(currUserHome, 0.7, 0.3);
    }

    private static void execution(String currUserHome, double alpha, double beta) throws IOException {
        System.out.printf("%.1f-%.1f\n", alpha, beta);
        String targetFileName = currUserHome +
                "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength_8_0.2_selection.bed";
        String actualResultFileName = String.format(
                "%s/experiments/ASM/simulation/CPGI_%.1f_%.1f/detection_summary_i90_r1_chr20_CPGI_%.1f_%.1f",
                currUserHome, alpha, beta, alpha, beta);
        String expectedResultFileName = String.format(
                "%s/experiments/ASM/simulation/CPGI_%.1f_%.1f/intersections_CPGI_%.1f_%.1f_4_5_10", currUserHome, alpha,
                beta, alpha, beta);

        String outputFileName = String.format(
                "%s/experiments/ASM/simulation/CPGI_%.1f_%.1f/evaluation_CPGI_%.1f_%.1f_4_5_10", currUserHome, alpha,
                beta, alpha, beta);

        twoWayIntersection(targetFileName, actualResultFileName, expectedResultFileName);
        fourWayIntersection(expectedResultFileName, actualResultFileName, outputFileName);
        System.out.println();
    }

    private static void twoWayIntersection(String targetFileName, String actualResultFileName,
                                           String expectedResultFileName) throws IOException {
        List<GenomicRegion> targetRegions = readBedRegions(targetFileName);
        List<GenomicRegion> resultRegions = readBedRegions(actualResultFileName);

        BufferedWriter writer = new BufferedWriter(new FileWriter(expectedResultFileName));

        for (GenomicRegion targetRegion : targetRegions) {
            for (GenomicRegion resultRegion : resultRegions) {
                if (targetRegion.getStart() <= resultRegion.getEnd() &&
                        targetRegion.getEnd() >= resultRegion.getStart()) {
                    targetRegion.setIntersected(true);
                    resultRegion.setIntersected(true);
                    writer.write(String.format("%s\t%s\t+\n", resultRegion.toBedString(), targetRegion.toBedString()));
                }
            }
        }

        for (GenomicRegion targetRegion : targetRegions) {
            if (!targetRegion.isIntersected()) {
                writer.write(targetRegion.toBedString() + "\t*" + "\n");
            }
        }

        for (GenomicRegion resultRegion : resultRegions) {
            if (!resultRegion.isIntersected()) {
                writer.write(resultRegion.toBedString() + "\t-" + "\n");
            }
        }

        writer.close();
    }

    private static void fourWayIntersection(String expectedResultFileName, String actualResultFileName,
                                            String outputFileName) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
        int tp = 0, fp = 0, fn = 0, tn = 0;
        List<GenomicRegion> expectedRegions = readBedRegions(expectedResultFileName, true);
        List<GenomicRegion> actualRegions = readBedRegions(actualResultFileName, true);
        for (GenomicRegion expectedRegion : expectedRegions) {
            for (GenomicRegion actualRegion : actualRegions) {
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
        System.out.println("accuracy:" + (tp + tn) / (double) (tp + tn + fn + fp));
        System.out.println("Precision:" + tp / (double) (tp + fp));
        System.out.println("Sensitivity:" + tp / (double) (tp + fn));

    }
}
