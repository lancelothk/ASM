package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

/**
 * Created by lancelothk on 11/14/14.
 * Utils for CPMR
 */
public class Utils {


    /**
     * write single interval output in ExtEpiRead format *
     */
    public static void writeExtEpireadInInterval(String outputPath, RefChr refChr, int startPos, int endPos,
                                                 Collection<MappedRead> mappedReadSet) throws IOException {
        BufferedWriter mappedReadWriter = new BufferedWriter(
                new FileWriter(String.format("%s/%s-%d-%d", outputPath, refChr.getChr(), startPos, endPos + 1)));
        for (MappedRead mappedRead : mappedReadSet) {
            if (mappedRead.getEpiRead(startPos, endPos) != null) {
                mappedReadWriter.write(mappedRead.getEpiRead(startPos, endPos).extEpireadFormat());
            }
        }
        mappedReadWriter.close();
    }

    /**
     * write single interval output in mapped read format like the original data *
     */
    public static void writeMappedReadInInterval(String outputPath, RefChr refChr, int startPos, int endPos,
                                                 Collection<MappedRead> mappedReadSet) throws IOException {
        BufferedWriter mappedReadWriter = new BufferedWriter(
                new FileWriter(String.format("%s/%s-%d-%d", outputPath, refChr.getChr(), startPos, endPos)));
        mappedReadWriter.write(String.format("ref:\t%s\n", refChr.getRefString().substring(startPos, endPos + 1)));
        for (MappedRead mappedRead : mappedReadSet) {
            mappedReadWriter.write(mappedRead.outputString(startPos, endPos));
        }
        mappedReadWriter.close();
    }
}
