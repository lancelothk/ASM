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
        List<GenomicRegion> targetRegions = readBedRegions(currUserHome +
                                                                   "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_qualifiedLength_8_0.2_selection.bed");
        List<GenomicRegion> resultRegions = readBedRegions(
                currUserHome + "/experiments/ASM/simulation/i90_r1_chr20_CPGI_1.0_0.0_detection_summary");

        BufferedWriter writer = new BufferedWriter(
                new FileWriter(currUserHome + "/experiments/ASM/simulation/intersections_CPGI_1_0_4_5_10"));

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
