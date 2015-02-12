package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.DataType.BedRegion;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import edu.cwru.cbc.ASM.commons.Utils;

import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;
import static edu.cwru.cbc.ASM.commons.Utils.readBedRegions;

/**
 * Created by kehu on 2/12/15.
 * generate semi-simulated ASM data
 */
public class Simulation {
    public static void main(String[] args) {

    }

    public static void executeSimulation(String referenceGenomeFileName, String readsFileName,
                                         String targetRegionFileName, double alpha, double beta) throws IOException {
        // read reference and refCpGs
        RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);

        // read target regions
        List<BedRegion> targetRegionList = readBedRegions(targetRegionFileName);


    }
}
