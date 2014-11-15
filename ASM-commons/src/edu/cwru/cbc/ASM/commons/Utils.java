package edu.cwru.cbc.ASM.commons;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;

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
//        return Files.asCharSource(intervalFile, Charsets.UTF_8).readLines(new LineProcessor<List<EpiRead>>() {
//            private List<EpiRead> epiReadList = new ArrayList<>();
//
//            @Override
//            public boolean processLine(String s) throws IOException {
//                epiReadList.add(new EpiRead(s, format));
//                return true;
//            }
//
//            @Override
//            public List<EpiRead> getResult() {
//                return epiReadList;
//            }
//        });
    }

}
