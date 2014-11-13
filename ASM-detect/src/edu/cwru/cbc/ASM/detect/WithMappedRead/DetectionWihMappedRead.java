package edu.cwru.cbc.ASM.detect.WithMappedRead;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MethylStatus;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;
import edu.cwru.cbc.ASM.detect.DataType.Edge;
import edu.cwru.cbc.ASM.detect.DataType.Vertex;
import edu.cwru.cbc.ASM.detect.Detection;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by kehu on 11/12/14.
 */
public class DetectionWihMappedRead extends Detection {
    private static final int MIN_READ_CPG = 2;
    private File intervalFile;
    private int initPos;
    private BufferedWriter groupWriter;


    public DetectionWihMappedRead(File intervalFile, int initPos, BufferedWriter groupWriter) {
        this.intervalFile = intervalFile;
        this.initPos = initPos;
        this.groupWriter = groupWriter;
    }

    public String call() throws Exception {
        // load input
        String reference = readRef(intervalFile);
        List<RefCpG> cpgList = extractCpGSite(reference, initPos);
        List<MappedRead> mappedReadList = Files.asCharSource(intervalFile, Charsets.UTF_8).readLines(
                new MappedReadFileLineProcessor());
        associateReadWithCpG(cpgList, mappedReadList);
        String[] items = intervalFile.getName().split("-");
        if (items.length != 3) {
            throw new RuntimeException("illegal interval file name format!");
        }

        StringBuilder asm_result = new StringBuilder(intervalFile.getName() + "\n");
        StringBuilder summary = new StringBuilder(intervalFile.getName());
        summary.append("\t").append(Integer.parseInt(items[2]) - Integer.parseInt(items[1]) + 1);
        summary.append("\t").append(mappedReadList.size());
        summary.append("\t").append(cpgList.size());

        List<Edge> edgeList = new ArrayList<>();
        Map<String, Vertex> vertexMap = new HashMap<>();

        constructGraph(vertexMap, edgeList, mappedReadList);
        System.out.println(
                intervalFile.getName() + "\t" + "vertex number:" + vertexMap.keySet().size() + "\t" + "edge number:" +
                        edgeList.size());
        getClusters(asm_result, edgeList, vertexMap);

        List<List<RefCpG>> groupCpGResults = new ArrayList<>();
        for (Vertex vertex : vertexMap.values()) {
            List<RefCpG> refCpGList = new ArrayList<>();
            asm_result.append("group: ").append(groupCpGResults.size()).append("\t").append(
                    vertex.getMappedReadList().size()).append("\n");
            vertex.getMappedReadList().sort(
                    (MappedRead read1, MappedRead read2) -> read1.getStart() - read2.getStart());
            for (MappedRead mappedRead : vertex.getMappedReadList()) {
                asm_result.append(mappedRead.getId()).append(",");
            }
            asm_result.append("\n");
            asm_result.append("Covered CpG sites:\n");
            List<CpG> vertexCpGList = new ArrayList<>(vertex.getCoveredCpGSites());
            vertexCpGList.sort((CpG c1, CpG c2) -> c1.getPos() - c2.getPos());
            for (CpG cpg : vertexCpGList) {
                RefCpG refCpG = new RefCpG(cpg.getPos());
                for (MappedRead mappedRead : vertex.getMappedReadList()) {
                    if (mappedRead.getFirstCpG().getPos() <= cpg.getPos() &&
                            mappedRead.getLastCpG().getPos() >= cpg.getPos()) {
                        if (mappedRead.getCpGMethylStatus(cpg.getPos()) == MethylStatus.C) {
                            refCpG.addMethylCount(1);
                        } else {
                            refCpG.addNonMethylCount(1);
                        }
                    }
                }
                asm_result.append(String.format("%d-%.2f,", cpg.getPos(), refCpG.getMethylLevel()));
                refCpGList.add(refCpG);
            }
            asm_result.append("\n");
            groupCpGResults.add(refCpGList);
        }
        asm_result.append("Number of groups:\t").append(vertexMap.values().size()).append("\n");
        summary.append("\t").append(vertexMap.values().size());
        groupWriter.write(asm_result.toString());


        Map<Integer, Integer> cpgPosMap = cpgList.stream().collect(Collectors.toMap(RefCpG::getPos, refCpG -> 0));
        for (Vertex vertex : vertexMap.values()) {
            for (CpG cpg : vertex.getCoveredCpGSites()) {
                if (cpgPosMap.containsKey(cpg.getPos())) {
                    cpgPosMap.put(cpg.getPos(), cpgPosMap.get(cpg.getPos()) + 1);
                }
            }
        }

        List<Integer> sortedPosList = new ArrayList<>(cpgPosMap.keySet());
        sortedPosList.sort(Integer::compareTo);
        double sum = 0;
        for (Integer integer : sortedPosList) {
            groupWriter.write(String.format("pos: %d\tcount: %d\n", integer, cpgPosMap.get(integer)));
            sum += cpgPosMap.get(integer);
        }

        Map<Integer, Integer> refCpGMap = new HashMap<>();
        for (int i = 0; i < sortedPosList.size(); i++) {
            // key is pos, value is index
            refCpGMap.put(sortedPosList.get(i), i);
        }

        writeAlignedGroupResult(refCpGMap, groupCpGResults, intervalFile);

        summary.append("\t").append(sum / sortedPosList.size()).append("\n");
        return summary.toString();
    }

    private void writeAlignedGroupResult(Map<Integer, Integer> refCpGMap, List<List<RefCpG>> groupCpGResult,
                                         File intervalFile) throws IOException {
        BufferedWriter alignedGroupWriter = new BufferedWriter(
                new FileWriter(intervalFile.getName() + ".alignedGroups"));
        // sort the cpg in each list first, then sort group list
        groupCpGResult.forEach(group -> group.sort(RefCpG::compareTo));
        groupCpGResult.sort((g1, g2) -> g1.get(0).getPos() - g2.get(0).getPos());

        int startIndex = refCpGMap.get(groupCpGResult.get(0).get(0).getPos());
        for (int i = 0; i < groupCpGResult.size(); i++) {
            alignedGroupWriter.write(i + "\t");
            int currIndex = refCpGMap.get(groupCpGResult.get(i).get(0).getPos());
            // fill the gap
            for (int j = 0; j < currIndex - startIndex; j++) {
                alignedGroupWriter.write(Strings.repeat(".", 20));
            }
            for (RefCpG refCpG : groupCpGResult.get(i)) {
                alignedGroupWriter.write(
                        Strings.padEnd(String.format("%f(%d)", refCpG.getMethylLevel(), refCpG.getCoveredCount()), 20,
                                       ' '));
            }
            alignedGroupWriter.write("\n");
        }
        alignedGroupWriter.close();
    }

    private void getClusters(StringBuilder asm_result, List<Edge> edgeList,
                             Map<String, Vertex> vertexMap) throws IOException {
        // count how many tie situation occurs. For analysis use
        int tieWeightCounter = 0, tieIdCountCounter = 0, tieMethylPolarity = 0;
        // merge vertexes connected by positive weight edge
        while (true) {
            // if edgeList is empty
            if (edgeList.size() == 0) {
                break;
            }
            List<Edge> maxEdgeList = getMaxFromList(edgeList, Edge::getWeight);
            // if max weight <= 0, stop merge.
            if (maxEdgeList.get(0).getWeight() < 0) {
                break;
            } else {
                // since edgeList is not empty, getMaxFromList always returns at least one element
                if (maxEdgeList.size() != 1) {
                    // multiple equal max weight edges.
                    tieWeightCounter++;
                    // use methyl polarity to group similar vertex first
                    // sort edge by abs of methyl polarity. Larger first
//                    maxEdgeList = getMaxFromList(maxEdgeList, Edge::getMethylPolarityAbs);
//                    if (maxEdgeList.size() != 1) {
//                        tieMethylPolarity++;

                    // First pick edge not connect to two clusters.
                    // sort edge by id count. Smaller first
                    maxEdgeList = getMinFromList(maxEdgeList, Edge::getIdCount);
                    if (maxEdgeList.size() != 1) {
                        tieIdCountCounter++;
                    }
//                    }
                }
                mergeVertex(maxEdgeList.get(0), edgeList, vertexMap);
            }
        }
        asm_result.append(String.format("tied weight counter:%d\n", tieWeightCounter));
        asm_result.append(String.format("tied methyl polarity:%d\n", tieMethylPolarity));
        asm_result.append(String.format("tied id counter:%d\n", tieIdCountCounter));
    }

    private void mergeVertex(Edge edge, List<Edge> edgeList, Map<String, Vertex> vertexMap) {
        Vertex left = edge.getLeft();
        Vertex right = edge.getRight();
        // merge right to left vertex
        left.addMappedRead(right.getMappedReadList());
        left.addCpG(right.getCoveredCpGSites());
        //remove this edge and right vertex from graph
        edge.removeFromVertex();
        edgeList.remove(edge);
        vertexMap.remove(right.getId());
        // update right vertex adj edges
        for (Edge edgeR : right.getAdjEdges()) {
            // update new vertex
            if (!edgeR.replaceVertex(right, left)) {
                throw new RuntimeException("fail to update new vertex in adj edge!");
            }
        }
        // update left vertex with new edges
        for (Edge adjEdge : right.getAdjEdges()) {
            left.addEdge(adjEdge);
        }
        // update edges of vertex which connect both left and right
        updateAndRemoveDupEdge(edgeList);
    }

    private void updateAndRemoveDupEdge(List<Edge> edgeList) {
        Map<String, Edge> edgeMap = new HashMap<>();
        Iterator<Edge> edgeIterator = edgeList.iterator();
        while (edgeIterator.hasNext()) {
            Edge edgeB = edgeIterator.next();
            // duplicate edge
            if (edgeMap.containsKey(edgeB.getUniqueId())) {
                Edge edgeA = edgeMap.get(edgeB.getUniqueId());
                // keep edgeA, remove edgeB.
                edgeA.setWeight(edgeA.getWeight() + edgeB.getWeight());
                edgeB.removeFromVertex();
                edgeIterator.remove(); // remove from edgelist
            } else {
                edgeMap.put(edgeB.getUniqueId(), edgeB);
            }
        }
    }

    /**
     * select items with max property value in the list.
     *
     * @return list contains items with max weight
     */
    private <T> List<T> getMaxFromList(List<T> itemList, Function<T, Integer> getProperty) {
        List<T> maxItemList = new ArrayList<>();
        int maxValue = Integer.MIN_VALUE;
        for (T item : itemList) {
            if (getProperty.apply(item) > maxValue) {
                maxItemList.clear();
                maxItemList.add(item);
                maxValue = getProperty.apply(item);
            } else if (getProperty.apply(item) == maxValue) {
                maxItemList.add(item);
            }
        }
        return maxItemList;
    }

    /**
     * select items with min property value in the list.
     *
     * @return list contains items with min weight
     */
    private <T> List<T> getMinFromList(List<T> itemList, Function<T, Integer> getProperty) {
        List<T> minItemList = new ArrayList<>();
        int minValue = Integer.MAX_VALUE;
        for (T item : itemList) {
            if (getProperty.apply(item) < minValue) {
                minItemList.clear();
                minItemList.add(item);
                minValue = getProperty.apply(item);
            } else if (getProperty.apply(item) == minValue) {
                minItemList.add(item);
            }
        }
        return minItemList;
    }

    private void constructGraph(Map<String, Vertex> vertexMap, List<Edge> edgeList, List<MappedRead> mappedReadList) {
        // visit all possible edges once.
        for (int i = 0; i < mappedReadList.size(); i++) {
            // only use reads which contains at least two Cpg sites
            if (mappedReadList.get(i).getCpgList().size() >= MIN_READ_CPG) {
                for (int j = i + 1; j < mappedReadList.size(); j++) {
                    int score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
                    if (score != Integer.MIN_VALUE) {
                        String idI = mappedReadList.get(i).getId(), idJ = mappedReadList.get(j).getId();
                        if (!vertexMap.containsKey(idI)) {
                            vertexMap.put(idI, new Vertex(mappedReadList.get(i)));
                        }
                        if (!vertexMap.containsKey(idJ)) {
                            vertexMap.put(idJ, new Vertex(mappedReadList.get(j)));
                        }
                        edgeList.add(new Edge(vertexMap.get(idI), vertexMap.get(idJ), score));
                    }
                }
            }
        }
    }

    private int checkCompatible(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            int score = 0;
            for (CpG cpgA : readA.getCpgList()) {
                for (CpG cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos()) {
                        if (cpgA.getMethylStatus() != cpgB.getMethylStatus()) {
                            score--;
                        } else {
                            score++;
                        }
                    }
                }
            }
            return score;
        } else {
            // don't have overlapped CpG, non-connected
            return Integer.MIN_VALUE;
        }
    }

    private String readRef(File intervalFile) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(intervalFile));
        String line = bufferedReader.readLine();
        String[] items = line.split("\t");
        if (items.length != 2) {
            throw new RuntimeException("invalid reference line in interval read file!");
        } else {
            return items[1].toUpperCase();
        }
    }

    private void associateReadWithCpG(List<RefCpG> cpgList, List<MappedRead> mappedReadList) {
        for (MappedRead read : mappedReadList) {
            for (RefCpG cpg : cpgList) {
                if (cpg.getPos() >= read.getStart() && cpg.getPos() + 1 <= read.getEnd()) {
                    read.addCpG(new CpG(cpg, read.getCpGMethylStatus(cpg.getPos())));
                }
            }
        }
    }

    private List<RefCpG> extractCpGSite(String reference, int initPos) {
        reference = reference.replace(" ", "");
        List<RefCpG> cpgList = new ArrayList<>();
        for (int i = 0; i < reference.length() - 1; i++) {
            if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
                cpgList.add(new RefCpG(initPos + i));
            }
        }
        return cpgList;
    }
}
