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
        double alpha = 0.8, beta = 0.2;
        List<GenomicRegion> targetRegions = readBedRegions(currUserHome +
                                                                   "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength_8_0.2_selection.bed");
        List<GenomicRegion> resultRegions = readBedRegions(
                String.format("%s/experiments/ASM/simulation/detection_summary_i90_r1_chr20_CPGI_%.1f_%.1f",
                              currUserHome, alpha, beta));

        BufferedWriter writer = new BufferedWriter(new FileWriter(
                String.format("%s/experiments/ASM/simulation/intersections_CPGI_%.1f_%.1f_4_5_10", currUserHome, alpha,
                              beta)));

        for (GenomicRegion targetRegion : targetRegions) {
            for (GenomicRegion resultRegion : resultRegions) {
                if (targetRegion.getStart() <= resultRegion.getEnd() &&
                        targetRegion.getEnd() >= resultRegion.getStart()) {
                    targetRegion.setIntersected(true);
                    resultRegion.setIntersected(true);
                    writer.write(String.format("*****\ttarget-%s\tresult-%s\n", targetRegion.toString(),
                                               resultRegion.toString()));
                }
            }
        }

        for (GenomicRegion targetRegion : targetRegions) {
            if (!targetRegion.isIntersected()) {
                writer.write("+++++\t" + targetRegion.toString() + "\n");
            }
        }

        for (GenomicRegion resultRegion : resultRegions) {
            if (!resultRegion.isIntersected()) {
                writer.write("-----\t" + resultRegion.toString() + "\n");
            }
        }

        writer.close();
    }
}
