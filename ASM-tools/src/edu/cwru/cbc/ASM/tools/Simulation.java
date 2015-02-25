package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.DataType.*;
import edu.cwru.cbc.ASM.commons.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;
import static edu.cwru.cbc.ASM.commons.Utils.readBedRegions;

/**
 * Created by kehu on 2/12/15.
 * generate semi-simulated ASM data
 */
public class Simulation {
    public static void main(String[] args) throws IOException {
      	String currUserHome = System.getProperty("user.home");
        executeSimulation(currUserHome + "/experiments/ASM/data/hg18_chr20.fa",
                          currUserHome + "/experiments/ASM/data/i90_r1_chr20",
                          currUserHome + "/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20.bed",
                          currUserHome + "/experiments/ASM/simulation/i90_r1_chr20_CPGI", 0.8, 0.05);
    }

    /**
     * execute data simulation
     *
     * @param alpha read will keep same pattern with chosen allele with probability alpha.
     * @param beta  methyl status of CpG will flip with probability beta
     * @throws IOException
     */
    public static void executeSimulation(String referenceGenomeFileName, String readsFileName,
                                         String targetRegionFileName, String outputFileName, double alpha,
                                         double beta) throws IOException {

		outputFileName = String.format("%s_%f_%f.sim", outputFileName, alpha, beta);
        // read reference and refCpGs
        RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);


        // read target regions
        List<GenomicRegion> targetRegions = readBedRegions(targetRegionFileName);
        Collections.sort(targetRegions);

        // generate non-ASM regions
        List<GenomicRegion> nonASMRegions = generateNonASMRegions(refChr, targetRegions);

        // attach refCpG to regions
        attachRefCpGToRegions(refCpGList, targetRegions, nonASMRegions);

        // read input sequences
        List<MappedRead> mappedReadList = Files.asCharSource(new File(readsFileName), Charsets.UTF_8).readLines(
                new MappedReadLineProcessor(refCpGList));

        // assign methyl status
        assignMethylStatusForNonASMRegion(nonASMRegions);
        assignMethylStatusForASMRegion(targetRegions);

        // add randomness
        addRandomnessToReads(nonASMRegions, alpha, beta);
        addRandomnessToReads(targetRegions, alpha, beta);

        // write out result
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
        for (MappedRead mappedRead : mappedReadList) {
            writer.write(mappedRead.toSimulationString() + "\n");
        }
        writer.close();
    }

    private static void addRandomnessToReads(List<GenomicRegion> regions, double alpha, double beta) {
        for (GenomicRegion region : regions) {
            Set<MappedRead> readsInRegion = getReadsInRegion(region);
            for (MappedRead read : readsInRegion) {
                Random rand = new Random();
                if (rand.nextDouble() > alpha) {
                    // add randomness
                    for (CpG cpG : read.getCpgList()) {
                        if (rand.nextDouble() < beta) {
                            // add randomness. Flip methyl status
                            cpG.setMethylStatus(
                                    cpG.getMethylStatus() == MethylStatus.C ? MethylStatus.T : MethylStatus.C);
                        }
                    }
                }
            }
        }
    }

    private static void assignMethylStatusForASMRegion(List<GenomicRegion> targetRegionList) {
        // each refCpG is attached to a region, each CpG in read belongs to a refCpG.
        // ASM region
        for (GenomicRegion targetRegion : targetRegionList) {
            Set<MappedRead> mappedReadSet = getReadsInRegion(targetRegion);

            // randomly assign read to each allele
            Random rand = new Random();
            for (MappedRead mappedRead : mappedReadSet) {
                if (rand.nextBoolean()) {
                    // allele 1
                    for (CpG cpg : mappedRead.getCpgList()) {
                        if (targetRegion.getRefMethylStatus(cpg.getRefCpG().getIndex())) {
                            cpg.setMethylStatus(MethylStatus.C);
                        } else {
                            cpg.setMethylStatus(MethylStatus.T);
                        }
                    }
                } else {
                    // allele 2 : conterpart of allele 1
                    for (CpG cpg : mappedRead.getCpgList()) {
                        if (targetRegion.getRefMethylStatus(cpg.getRefCpG().getIndex())) {
                            cpg.setMethylStatus(MethylStatus.T); // opposite to allele 1
                        } else {
                            cpg.setMethylStatus(MethylStatus.C); // opposite to allele 1
                        }
                    }
                }
            }
        }
    }

    private static Set<MappedRead> getReadsInRegion(GenomicRegion region) {
        Set<MappedRead> mappedReadSet = new HashSet<>();
        for (int i = 0; i < region.getRefCpGList().size(); i++) {
            RefCpG refCpG = region.getRefCpGList().get(i);
            refCpG.assignIndex(i);
            for (CpG cpg : refCpG.getCpGList()) {
                mappedReadSet.add(cpg.getMappedRead());
            }
        }
        return mappedReadSet;
    }

    private static void assignMethylStatusForNonASMRegion(List<GenomicRegion> nonASMRegions) {
        // each refCpG is attached to a region, each CpG in read belongs to a refCpG.
        // Non ASM region
        for (GenomicRegion nonASMRegion : nonASMRegions) {
            for (int i = 0; i < nonASMRegion.getRefCpGList().size(); i++) {
                for (CpG cpg : nonASMRegion.getRefCpGList().get(i).getCpGList()) {
                    if (nonASMRegion.getRefMethylStatus(i)) {
                        cpg.setMethylStatus(MethylStatus.C);
                    } else {
                        cpg.setMethylStatus(MethylStatus.T);
                    }
                }
            }
        }
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
        for (GenomicRegion region : allRegions) {
            region.generateRandomAllelePattern();
        }
    }

    /**
     * generate Non ASM regions which covers the region not covered by ASM regions.
     *
     * @return list of non ASM regions. In position increasing order.
     */
    private static List<GenomicRegion> generateNonASMRegions(RefChr refChr, List<GenomicRegion> targetRegionList) {
        List<GenomicRegion> nonASMRegions = new ArrayList<>();
        int lastEnd = -1;
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
