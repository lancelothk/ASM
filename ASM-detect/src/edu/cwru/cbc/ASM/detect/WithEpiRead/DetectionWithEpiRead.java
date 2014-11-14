package edu.cwru.cbc.ASM.detect.WithEpiRead;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;
import edu.cwru.cbc.ASM.commons.Utils;
import edu.cwru.cbc.ASM.detect.Detection;
import edu.cwru.cbc.ASM.detect.WithEpiRead.DataType.Edge;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by kehu on 11/12/14.
 * Implementation of ASM detection with using EpiRead
 */
public class DetectionWithEpiRead extends Detection {
    private static final int MIN_READ_CPG = 2;
    private File intervalFile;
    private BufferedWriter groupWriter;
    private EpiRead.EpiReadFormat format;

    public DetectionWithEpiRead(File intervalFile, BufferedWriter summaryWriter, EpiRead.EpiReadFormat format) {
        this.intervalFile = intervalFile;
        this.groupWriter = summaryWriter;
        this.format = format;
    }

    private double[][] buildInputMatrix() throws IOException {
        List<EpiRead> epireadList = Utils.readEpiReadFile(intervalFile, format);
        Map<Integer, Integer> cpgMap = new HashMap<>();
        epireadList.forEach(epiRead -> {
            if (!cpgMap.containsKey(epiRead.getCpgOrder())) {
                cpgMap.put(epiRead.getCpgOrder(), epiRead.getCpgPos());
            }
        });
        List<Integer> cpgs = new ArrayList<>(cpgMap.keySet());
        cpgs.sort(Integer::compare);
        double[][] inputMatrix = new double[epireadList.size()][cpgs.size()];
        for (int i = 0; i < inputMatrix.length; i++) {
            EpiRead epiread = epireadList.get(i);
            int firstIndex = cpgs.indexOf(epiread.getCpgOrder());
            for (int j = 0; j < epiread.getCpgSeq().length(); j++) {
                if (epiread.getCpgSeq().charAt(j) == 'C') {
                    inputMatrix[i][j + firstIndex] = 1;
                } else if (epiread.getCpgSeq().charAt(j) == 'T') {
                    inputMatrix[i][j + firstIndex] = -1;
                }
            }
        }
        return inputMatrix;
    }

//    private double[][] buildConflictMatrix(double[][] inputMatrix){
//        double[][] conflictMatrix = new double[inputMatrix.length][inputMatrix.length];
//        for (int i = 0; i < inputMatrix.length; i++) {
//            for (int j = i+1; j < inputMatrix.length; j++) {
//                for (int k = 0; k < inputMatrix[0].length; k++) {
//                    if (inputMatrix[i][k] == inputMatrix[j][k]){
//                        conflictMatrix[i][j] = 0;
//                    } else if (inputMatrix[i][k] != inputMatrix[j][k] && inputMatrix[i][k] != 0 && inputMatrix[j][k] != 0){
//                        conflictMatrix[i][j] = 1;
//                    }else {
//                        conflictMatrix[i][j] = 0.5;
//                    }
//                }
//            }
//        }
//        return conflictMatrix;
//    }

    @Override
    public String call() throws Exception {
        double[][] inputMatrix = buildInputMatrix();
//        double[][] conflictMatrix = buildConflictMatrix(inputMatrix);
        List<Edge> edgeList = new ArrayList<>();
        for (int i = 0; i < inputMatrix.length; i++) {
            for (int j = i + 1; j < inputMatrix.length; j++) {
                int score = 0, coveredCount = 0;
                for (int k = 0; k < inputMatrix[0].length; k++) {
                    if (inputMatrix[i][k] == inputMatrix[j][k]) {

                    } else if (inputMatrix[i][k] != inputMatrix[j][k] && inputMatrix[i][k] != 0 &&
                            inputMatrix[j][k] != 0) {
                        edgeList.add(new Edge(i, j, 1));
                    } else {
                        conflictMatrix[i][j] = 0.5;
                    }
                }
                edgeList.add(new Edge(i, j, 0));
            }
        }
        return "";
    }
}
