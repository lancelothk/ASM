package edu.cwru.cbc.ASM.detect.WithMappedRead.DataType;

import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MethylStatus;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;

/**
 * Created by kehu on 11/17/14.
 * Graph for detecting ASM.
 */
public class ASMGraph {
    private List<Edge> edgeList = new ArrayList<>();
    private Map<String, Vertex> vertexMap = new HashMap<>();
    private Map<Integer, Integer> coveredCpGMap; // pos <--> count
    private Map<String, Vertex> clusterResult;
    // count how many tie situation occurs. For analysis use
    private int tieWeightCounter, tieIdCountCounter, tieMethylPolarity;
    private int originalVertexCount, originalEdgeCount;

    public ASMGraph(List<MappedRead> mappedReadList) {
        // visit all possible edges once.
        for (int i = 0; i < mappedReadList.size(); i++) {
            for (int j = i + 1; j < mappedReadList.size(); j++) {
                double score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
                if (score != Double.MIN_VALUE) {
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
        this.originalVertexCount = vertexMap.size();
        this.originalEdgeCount = edgeList.size();
    }

    public void getClusters() throws IOException {
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
                mergeVertex(maxEdgeList.get(0));
            }
        }
        this.clusterResult = vertexMap;
        setCoveredCpGMap();
    }

    private void mergeVertex(Edge edge) {
        Vertex left = edge.getLeft();
        Vertex right = edge.getRight();
        // merge right to left vertex
        left.addMappedRead(right.getMappedReadList());
        left.addCpG(right.getRefCpGMap().values());
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
        right.getAdjEdges().forEach(left::addEdge);
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
    private <T> List<T> getMaxFromList(List<T> itemList, Function<T, Double> getProperty) {
        List<T> maxItemList = new ArrayList<>();
        double maxValue = getProperty.apply(itemList.get(0));
        for (T item : itemList) {
            if (getProperty.apply(item) > maxValue) {
                maxItemList.clear();
                maxItemList.add(item);
                maxValue = getProperty.apply(item);
            } else if (getProperty.apply(item) == maxValue) {
                maxItemList.add(item);
            }
        }
        if (maxItemList.size() == 0) {
            throw new RuntimeException("max item list is empty!");
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
        int minValue = getProperty.apply(itemList.get(0));
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

    private double checkCompatible(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            int score = 0, count = 0;
            for (CpG cpgA : readA.getCpgList()) {
                for (CpG cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos()) {
                        if (cpgA.getMethylStatus() != cpgB.getMethylStatus()) {
                            if (cpgA.getMethylStatus() == MethylStatus.N || cpgB.getMethylStatus() == MethylStatus.N) {
                                score -= 0.5;
                            } else {
                                score--;
                            }
                        } else {
                            score++;
                        }
                        count++;
                    }
                }
            }
            return score / (double) count;
        } else {
            // don't have overlapped CpG, non-connected
            return 0;
        }
    }

    private void setCoveredCpGMap() {
        coveredCpGMap = new HashMap<>();
        for (Vertex vertex : vertexMap.values()) {
            for (RefCpG refCpG : vertex.getRefCpGMap().values()) {
                if (!coveredCpGMap.containsKey(refCpG.getPos())) {
                    coveredCpGMap.put(refCpG.getPos(), 1);
                } else {
                    coveredCpGMap.put(refCpG.getPos(), coveredCpGMap.get(refCpG.getPos()) + 1);
                }
            }
        }
    }

    public Map<Integer, Integer> getCoveredCpGMap() {
        return coveredCpGMap;
    }

    public Map<String, Vertex> getClusterResult() {
        return clusterResult;
    }

    public int getTieWeightCounter() {
        return tieWeightCounter;
    }

    public int getTieIdCountCounter() {
        return tieIdCountCounter;
    }

    public int getTieMethylPolarity() {
        return tieMethylPolarity;
    }

    public int getOriginalVertexCount() {
        return originalVertexCount;
    }

    public int getOriginalEdgeCount() {
        return originalEdgeCount;
    }
}
