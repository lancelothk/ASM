package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.DataType.GenomicRegion;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import edu.cwru.cbc.ASM.commons.Utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;
import static edu.cwru.cbc.ASM.commons.Utils.readBedRegions;

/**
 * Created by kehu on 2/12/15.
 * generate semi-simulated ASM data
 */
public class Simulation {
    public static void main(String[] args) throws IOException {
        executeSimulation("/home/kehu/experiments/ASM/data/hg18_chr20.fa", "",
                          "/home/kehu/experiments/ASM/simulation/simu.bed", 0, 0);
    }

    public static void executeSimulation(String referenceGenomeFileName, String readsFileName,
                                         String targetRegionFileName, double alpha, double beta) throws IOException {
        // read reference and refCpGs
        RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);


        // read target regions
        List<GenomicRegion> targetRegionList = readBedRegions(targetRegionFileName);
        targetRegionList.sort(new Comparator<GenomicRegion>() {
            @Override
            public int compare(GenomicRegion o1, GenomicRegion o2) {
                int startCompare = Long.compare(o1.getStart(), o2.getStart());
                int endCompare = Long.compare(o1.getEnd(), o2.getEnd());
                return startCompare == 0 ? endCompare : startCompare;
            }
        });

        // generate non-ASM regions
        List<GenomicRegion> nonASMRegions = generateNonASMRegions(refChr, targetRegionList);


        // attach refCpG to regions
        attachRefCpGToRegions(refCpGList, targetRegionList, nonASMRegions);


    }

    private static void attachRefCpGToRegions(List<RefCpG> refCpGList, List<GenomicRegion> targetRegionList,
                                              List<GenomicRegion> nonASMRegions) {
        List<GenomicRegion> allRegions = new ArrayList<>();
        allRegions.addAll(nonASMRegions);
        allRegions.addAll(targetRegionList);
        for (RefCpG refCpG : refCpGList) {
            for (GenomicRegion region : allRegions) {
                if (refCpG.getPos() >= region.getStart() && refCpG.getPos() <= region.getEnd()) {
                    region.addRefCpG(refCpG);
                    break;
                }
            }
        }
    }

    /**
     * generate Non ASM regions which covers the region not covered by ASM regions.
     *
     * @return list of non ASM regions. In position increasing order.
     */
    private static List<GenomicRegion> generateNonASMRegions(RefChr refChr, List<GenomicRegion> targetRegionList) {
        List<GenomicRegion> nonASMRegions = new ArrayList<>();
        long lastEnd = -1;
        for (int i = 0; i <= targetRegionList.size(); i++) {
            if (i == 0) {
                nonASMRegions.add(
                        new GenomicRegion(refChr.getChr(), refChr.getStart(), targetRegionList.get(i).getStart() - 1,
                                          String.valueOf(nonASMRegions.size() + 1)));
                lastEnd = targetRegionList.get(i).getEnd();
            } else if (i == targetRegionList.size()) {
                nonASMRegions.add(new GenomicRegion(refChr.getChr(), lastEnd + 1, refChr.getEnd(),
                                                    String.valueOf(nonASMRegions.size() + 1)));
            } else {
                nonASMRegions.add(
                        new GenomicRegion(refChr.getChr(), lastEnd + 1, targetRegionList.get(i).getStart() - 1,
                                          String.valueOf(nonASMRegions.size() + 1)));
                lastEnd = targetRegionList.get(i).getEnd();
            }
        }
        return nonASMRegions;
    }
}
