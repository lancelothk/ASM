package edu.cwru.cbc.ASM.commons;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 11/13/14.
 */
public class Utils {

    public static List<EpiRead> readEpiReadFile(File inputFile, EpiRead.EpiReadFormat format) throws IOException {
        List<EpiRead> epiReadList = new ArrayList<>();
        Files.lines(inputFile.toPath()).forEach(line -> epiReadList.add(EpiRead.ParseEpiRead(line, format)));
        return epiReadList;
    }

}
