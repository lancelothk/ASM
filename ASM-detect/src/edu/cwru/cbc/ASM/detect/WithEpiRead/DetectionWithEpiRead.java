package edu.cwru.cbc.ASM.detect.WithEpiRead;

import edu.cwru.cbc.ASM.commons.DataType.EpiRead;
import edu.cwru.cbc.ASM.commons.Utils;
import edu.cwru.cbc.ASM.detect.Detection;
import edu.cwru.cbc.ASM.detect.WithEpiRead.DataType.Edge;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.logging.Logger;

/**
 * Created by kehu on 11/12/14.
 * Implementation of ASM detection with using EpiRead
 */
public class DetectionWithEpiRead extends Detection {
    private static final Logger logger = Logger.getLogger(DetectionWithEpiRead.class.getName());
    private EpiRead.EpiReadFormat format;

    public DetectionWithEpiRead(EpiRead.EpiReadFormat format) {
        this.format = format;
    }

    public static void main(String[] args) throws ExecutionException, InterruptedException {
        long start = System.currentTimeMillis();
//        String pathName = "/home/lancelothk/experiment/ASM/epiread/intervals_epiread/chr22-14476217-14476662";
        String inputName = "/home/lancelothk/IdeaProjects/ASM/ASM-detect/testData/epiTest";
        int threadNumber = 6;
        EpiRead.EpiReadFormat format = EpiRead.EpiReadFormat.extEpiread;

        new DetectionWithEpiRead(format).execute(inputName, threadNumber);
        logger.finest(System.currentTimeMillis() - start + "ms");
    }

    private double[][] buildInputMatrix() throws IOException {
        List<EpiRead> epireadList = Utils.readEpiReadFile(inputFile, format);
        Map<Integer, Integer> cpgMap = new HashMap<>();
        epireadList.forEach(epiRead -> {
            int order = epiRead.getCpgOrder();
            for (int i = 0; i < epiRead.getCpgSeq().length(); i++) {
                if (!cpgMap.containsKey(order + i)) {
                    cpgMap.put(order + i, 0);
                }
            }
        });
        List<Integer> cpgs = new ArrayList<>(cpgMap.keySet());
        cpgs.sort(Integer::compareTo);
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
        extractIntervalPosition(inputFile);

        Map<Integer, Integer> p = new HashMap<>();
        Map<Integer, Integer> partA = new HashMap<>();
        Map<Integer, Integer> partB = new HashMap<>();
        double[][] inputMatrix = buildInputMatrix();
//        double[][] conflictMatrix = buildConflictMatrix(inputMatrix);
        List<Edge> edgeList = new ArrayList<>();
        for (int i = 0; i < inputMatrix.length; i++) {
            for (int j = i + 1; j < inputMatrix.length; j++) {
                double score = 0, coveredCount = 0;
                for (int k = 0; k < inputMatrix[0].length; k++) {
                    if (inputMatrix[i][k] != 0 || inputMatrix[j][k] != 0) {
                        coveredCount++;
                    }
                    if (inputMatrix[i][k] != inputMatrix[j][k]) {
                        if (inputMatrix[i][k] != 0 && inputMatrix[j][k] != 0) {
                            score++;
                        } else {
                            score += 0.5;
                        }
                    }
                }
                edgeList.add(new Edge(i, j, score / coveredCount));
            }
        }
        edgeList.sort(Edge::compareTo);
        int tieCount = 0;
        if (edgeList.size() >= 2 && edgeList.get(0).getWeight() == edgeList.get(1).getWeight()) {
            tieCount++;
        }
        // initialize partition
        partA.put(edgeList.get(0).getLeft(), 0);
        partB.put(edgeList.get(0).getRight(), 0);
        p.put(edgeList.get(0).getLeft(), 0);
        p.put(edgeList.get(0).getRight(), 0);
        edgeList.remove(0);

        while (edgeList.size() != 0) {
            for (int i = 0; i < edgeList.size(); i++) {
                Edge e = edgeList.get(i);
                if (partA.containsKey(e.getLeft()) || partA.containsKey(e.getRight()) ||
                        partB.containsKey(e.getLeft()) || partB.containsKey(e.getRight())) {
                    if (edgeList.size() >= 2 && edgeList.get(0).getWeight() == edgeList.get(1).getWeight()) {
                        tieCount++;
                    }
                    if (partA.containsKey(e.getLeft()) && !p.containsKey(e.getRight())) {
                        partB.put(e.getRight(), 0);
                        p.put(e.getRight(), 0);
                    } else if (partA.containsKey(e.getRight()) && !p.containsKey(e.getLeft())) {
                        partB.put(e.getLeft(), 0);
                        p.put(e.getLeft(), 0);
                    } else if (partB.containsKey(e.getLeft()) && !p.containsKey(e.getRight())) {
                        partA.put(e.getRight(), 0);
                        p.put(e.getRight(), 0);
                    } else if (partB.containsKey(e.getRight()) && !p.containsKey(e.getLeft())) {
                        partA.put(e.getLeft(), 0);
                        p.put(e.getLeft(), 0);
                    }
                    edgeList.remove(i);
                    break;
                }
            }
            for (int i = edgeList.size() - 1; i > 0; i--) {
                Edge e = edgeList.get(i);
                if (partA.containsKey(e.getLeft()) || partA.containsKey(e.getRight()) ||
                        partB.containsKey(e.getLeft()) || partB.containsKey(e.getRight())) {
                    if (edgeList.size() >= 2 && edgeList.get(0).getWeight() == edgeList.get(1).getWeight()) {
                        tieCount++;
                    }
                    if (partA.containsKey(e.getLeft()) && !p.containsKey(e.getRight())) {
                        partA.put(e.getRight(), 0);
                        p.put(e.getRight(), 0);
                    } else if (partA.containsKey(e.getRight()) && !p.containsKey(e.getLeft())) {
                        partA.put(e.getLeft(), 0);
                        p.put(e.getLeft(), 0);
                    } else if (partB.containsKey(e.getLeft()) && !p.containsKey(e.getRight())) {
                        partB.put(e.getRight(), 0);
                        p.put(e.getRight(), 0);
                    } else if (partB.containsKey(e.getRight()) && !p.containsKey(e.getLeft())) {
                        partB.put(e.getLeft(), 0);
                        p.put(e.getLeft(), 0);
                    }
                    edgeList.remove(i);
                    break;
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("PartA\n");
        partA.forEach((i, j) -> sb.append(i).append("\t"));
        sb.append("\n");

        sb.append("PartB\n");
        partB.forEach((i, j) -> sb.append(i).append("\t"));
        sb.append("\n");
        logger.info(sb.toString());
        return "";
    }

    @Override
    protected Detection constructNewInstance() {
        return new DetectionWithEpiRead(format);
    }
}
