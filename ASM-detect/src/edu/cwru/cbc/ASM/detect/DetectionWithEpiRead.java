package edu.cwru.cbc.ASM.detect;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by kehu on 11/12/14.
 */
public class DetectionWithEpiRead extends Detection {
    private static final int MIN_READ_CPG = 2;
    private File intervalFile;
    private BufferedWriter groupWriter;
    private EpiRead.EpiReadFormat format;

    public DetectionWithEpiRead(File intervalFile, BufferedWriter summaryWriter) {
        this.intervalFile = intervalFile;
        this.groupWriter = summaryWriter;
    }

    private List<EpiRead> loadInterval() throws IOException {
        List<EpiRead> epiReadList = new ArrayList<>();
        Files.lines(intervalFile.toPath()).forEach(line -> epiReadList.add(new EpiRead(line, format)));
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

    private void buildMatrix() throws IOException {
        List<EpiRead> epireadList = loadInterval();
        Map<Integer, Integer> cpgMap = new HashMap<>();
        epireadList.forEach(epiRead -> {
            if (!cpgMap.containsKey(epiRead.getCpgOrder())) {
                cpgMap.put(epiRead.getCpgOrder(), epiRead.getCpgPos());
            }
        });
        List<Integer> cpgs = new ArrayList<>(cpgMap.keySet());
        cpgs.sort(Integer::compare);
        int[][] matrix = new int[epireadList.size()][cpgs.size()];
        for (int i = 0; i < matrix.length; i++) {
            EpiRead epiread = epireadList.get(i);
            int firstIndex = cpgs.indexOf(epiread.getCpgOrder());
            for (int j = 0; j < epiread.getCpgSeq().length(); j++) {
                if (epiread.getCpgSeq().charAt(j) == 'C') {
                    matrix[i][j + firstIndex] = 1;
                } else if (epiread.getCpgSeq().charAt(j) == 'T') {
                    matrix[i][j + firstIndex] = -1;
                }
            }
        }
    }

    @Override
    public String call() throws Exception {
        buildMatrix();
        return "";
    }
}
