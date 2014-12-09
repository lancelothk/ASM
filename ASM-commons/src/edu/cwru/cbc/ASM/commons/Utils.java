package edu.cwru.cbc.ASM.commons;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 11/13/14.
 * Utils for methods shared by modules under ASM project.
 */
public class Utils {

    public static List<EpiRead> readEpiReadFile(File inputFile, EpiRead.EpiReadFormat format) throws IOException {
        List<EpiRead> epiReadList = new ArrayList<>();
        Files.lines(inputFile.toPath()).forEach(line -> epiReadList.add(EpiRead.ParseEpiRead(line, format)));
        return epiReadList;
    }

    /**
     * Extract reference CpG from reference string.
     * @param reference should be upperCase string and without space in it.
     * @return returned List is sorted by position.
     */
    public static List<RefCpG> extractCpGSite(String reference, int initPos) {
        // initialize refCpGList with 1/20 of reference size, since the probability of CG occurrence should be 1/16.
        // Actual probability is higher than 1/16 from current observation. In CPG dense region, it should be much higher.
        List<RefCpG> reFCpGList = new ArrayList<>(reference.length() / 20);
        for (int i = 0; i < reference.length() - 1; i++) {
            if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
                reFCpGList.add(new RefCpG(i + initPos));
            }
        }
        return reFCpGList;
    }

}
