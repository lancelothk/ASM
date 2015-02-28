package edu.cwru.cbc.ASM.CPMR;

import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.RefChr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

/**
 * Created by lancelothk on 11/14/14.
 * CommonsUtils for CPMR
 */
public class CPMRUtils {

    /**
     * write single interval output in mapped read format like the original data *
     */
    public static void writeMappedReadInInterval(String intervalFolderName, RefChr refChr, int startPos, int endPos,
                                                 Collection<MappedRead> mappedReadSet) throws IOException {
        BufferedWriter mappedReadWriter = new BufferedWriter(
                new FileWriter(String.format("%s/%s-%d-%d", intervalFolderName, refChr.getChr(), startPos, endPos)));
        mappedReadWriter.write(String.format("ref:\t%s\n", refChr.getRefString().substring(startPos, endPos + 1)));
        for (MappedRead mappedRead : mappedReadSet) {
            mappedReadWriter.write(mappedRead.outputString(startPos, endPos));
        }
        mappedReadWriter.close();
    }
}
