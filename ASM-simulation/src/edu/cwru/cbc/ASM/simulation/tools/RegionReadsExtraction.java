package edu.cwru.cbc.ASM.simulation.tools;

import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.GenomicInterval;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;
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

    public static final int MIN_READ_CPG = 2;

    public static void main(String[] args) throws IOException {
        long startTime = System.currentTimeMillis();
        String currUserHome = System.getProperty("user.home");
        String referenceGenomeFileName = currUserHome + "/experiments/ASM/data/hg18_chr20.fa";
        String inputReadsPath = currUserHome + "/experiments/ASM/data/i90_r1_chr20";
        String regionType = "CGI";
        String targetRegionFile = currUserHome +
                "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_" + regionType + ".bed";
        String outputFolder = currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/" + regionType + "_regions";

        File outputFolderFile = new File(outputFolder);
        if (!outputFolderFile.exists()) {
            outputFolderFile.mkdir();
        }

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
            String[] items = r.split("\t");
            int start = Integer.parseInt(items[2]);
            int end = Integer.parseInt(items[3]);

            MappedRead mappedRead = new MappedRead(items[0], items[1].charAt(0), start, end, items[4], items[5]);
            if (mappedRead.countCpG(refMap) >= MIN_READ_CPG) {
                targetRegionsMap.keySet()
                        .stream()
                        .filter(genomicInterval -> start <= genomicInterval.getEnd() &&
                                end >= genomicInterval.getStart())
                        .forEach(genomicInterval -> targetRegionsMap.get(genomicInterval).add(r));
            }

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
