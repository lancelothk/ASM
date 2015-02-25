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
        List<GenomicRegion> targetRegions = readBedRegions(
                "/home/kehu/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20.bed");
        List<GenomicRegion> resultRegions = readBedRegions(
                "/home/kehu/experiments/ASM/simulation/i90_r1_chr20_sim_CPGI_1_0_detection_summary");

        BufferedWriter writer = new BufferedWriter(
                new FileWriter("/home/kehu/experiments/ASM/simulation/intersections_CPGI_1_0_4_5_10_new"));

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
