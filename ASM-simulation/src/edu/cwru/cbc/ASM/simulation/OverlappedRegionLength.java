package edu.cwru.cbc.ASM.simulation;

import java.io.IOException;

/**
 * Created by kehu on 3/23/15.
 * Used to calculate overlapped region length
 */
public class OverlappedRegionLength {
    public static void main(String[] args) throws IOException {
        //        IntersectRegions.calcOverlappedLenth(
        //                "/home/kehu/experiments/ASM/simulation/CpGIslandsRegions/CpGIslandsRegions_selected_15_5.bed",
        //                "/home/kehu/experiments/ASM/simulation/test/0.7_0.3/experiment_4_2_5_10/result_summary_", "fn");
        IntersectRegions.calcOverlappedLenth(
                "/home/kehu/experiments/ASM/simulation/CpGIslandsRegions/CpGIslandsRegions_selected_15_5.bed",
                "/home/kehu/experiments/ASM/amrfinder/experiments/LRT_g1_w5_C0.5.amr");
    }
}
