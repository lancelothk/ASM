package edu.cwru.cbc.ASM.commons;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.BedRegion;
import edu.cwru.cbc.ASM.commons.DataType.EpiRead;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 11/13/14.
 * Utils for methods shared by modules under ASM project.
 */
public class Utils {

    public static List<EpiRead> readEpiReadFile(File inputFile, EpiRead.EpiReadFormat format) throws IOException {
        List<EpiRead> epiReadList = new ArrayList<>();
        java.nio.file.Files.lines(inputFile.toPath()).forEach(
                line -> epiReadList.add(EpiRead.ParseEpiRead(line, format)));
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

    public static RefChr readReferenceGenome(String inputFileName) throws IOException {
        List<String> lines = Files.readLines(new File(inputFileName), Charsets.UTF_8);
        StringBuilder referenceBuilder = new StringBuilder();
        // parse and remove first line, which is chromosome name
        String chr = lines.get(0).replace(">", "");
        lines.remove(0);
        lines.forEach(line -> referenceBuilder.append(line.replaceAll(" ", "").toUpperCase()));
        return new RefChr(chr, referenceBuilder.toString());
    }

    /**
     * read bed format regions
     *
     * @param bedFileName name of the input file
     * @return list of bed regions
     * @throws IOException
     */
    public static List<BedRegion> readBedRegions(String bedFileName) throws IOException {
        if (!bedFileName.endsWith(".bed")) {
            throw new RuntimeException("the input should be bed format file with .bed extension!");
        }
        return Files.readLines(new File(bedFileName), Charsets.UTF_8, new LineProcessor<List<BedRegion>>() {
            private List<BedRegion> bedRegionList = new ArrayList<>();

            @Override
            public boolean processLine(String line) throws IOException {
                String[] items = line.split("\t");
                if (items.length < 4) {
                    throw new RuntimeException("invalid bed format:" + line);
                }
                bedRegionList.add(
                        new BedRegion(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]), items[3]));
                return true;
            }

            @Override
            public List<BedRegion> getResult() {
                return bedRegionList;
            }
        });
    }
}
