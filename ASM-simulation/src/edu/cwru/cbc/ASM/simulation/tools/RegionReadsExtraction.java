package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.GenomicInterval;
import edu.cwru.cbc.ASM.commons.bed.BedUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/16/15.
 * extract reads by regions.
 */
public class RegionReadsExtraction {
    public static void main(String[] args) throws IOException {
        long startTime = System.currentTimeMillis();
        String currUserHome = System.getProperty("user.home");
        String referenceGenomeFileName = currUserHome + "/experiments/ASM/data/hg18_chr20.fa";
        String inputReadsPath = currUserHome + "/experiments/ASM/data/i90_r1_chr20";
        String targetRegionFile = currUserHome +
                "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_nonCGI.bed";
        String outputFolder = currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/nonCGI_regions";

        BufferedReader readReader = new BufferedReader(new FileReader(inputReadsPath));
        List<String> lines = readReader.lines().collect(Collectors.toList());
        readReader.close();
        System.out.println("finished loading reads");
        System.out.println(System.currentTimeMillis() - startTime + "ms");

        RefChr refChr = CommonsUtils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = CommonsUtils.extractCpGSite(refChr.getRefString(), 0);
        Map<Integer, RefCpG> refMap = refCpGList.stream().collect(Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));

        List<GenomicInterval> targetRegions = BedUtils.readSingleChromBedRegions(
                targetRegionFile);
        Map<GenomicInterval, List<String>> targetRegionsMap = new HashMap<>();
        for (GenomicInterval targetRegion : targetRegions) {
            targetRegionsMap.put(targetRegion, new ArrayList<>());
        }

        lines.stream().forEach(r -> {
            String[] tmpItems = r.split("\t");
            int x = Integer.parseInt(tmpItems[2]);
            int y = Integer.parseInt(tmpItems[3]);

            targetRegionsMap.keySet()
                    .stream()
                    .filter(genomicInterval -> x <= genomicInterval.getEnd() && y >= genomicInterval.getStart())
                    .forEach(genomicInterval -> targetRegionsMap.get(genomicInterval).add(r));

        });

        for (GenomicInterval genomicInterval : targetRegionsMap.keySet()) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(
                    String.format(outputFolder + "/i90_r1_chr20_" +
                            "%d-%d", genomicInterval.getStart(), genomicInterval.getEnd())));
            for (String r : targetRegionsMap.get(genomicInterval)) {
                writer.write(r + "\n");
            }
            writer.close();
        }

        System.out.println(System.currentTimeMillis() - startTime + "ms");
    }
}
