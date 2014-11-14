package edu.cwru.cbc.ASM.detect.WithMappedRead.DataType;

import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;

import java.util.*;

/**
 * Created by lancelothk on 5/14/14.
 * Vertex in graph.
 */
public class Vertex {
    private String id;
    private List<Edge> adjEdges;
    private List<MappedRead> mappedReadList;
    private Map<Integer, CpG> cpGMap;
    private int methylPolarity;

    public Vertex(MappedRead mappedRead) {
        this.id = mappedRead.getId();
        this.adjEdges = new ArrayList<>();
        this.mappedReadList = new ArrayList<>();
        this.mappedReadList.add(mappedRead);
        this.cpGMap = new HashMap<>();
        this.methylPolarity = mappedRead.getMethylPolarity();
        addCpG(mappedRead.getCpgList());
    }

    public void addCpG(Collection<CpG> cpgList) {
        for (CpG cpg : cpgList) {
            if (!cpGMap.containsKey(cpg.getPos())) {
                cpGMap.put(cpg.getPos(), cpg);
            }
        }
    }

    public Collection<CpG> getCoveredCpGSites() {
        return cpGMap.values();
    }

    public int getFirstCpGPos() {
        return cpGMap.values().stream().min((c1, c2) -> c1.getPos() - c2.getPos()).get().getPos();
    }

    public void addEdge(Edge edge) {
        this.adjEdges.add(edge);
    }

    public void removeEdge(Edge edge) {
        this.adjEdges.remove(edge);
    }

    public void addMappedRead(List<MappedRead> mappedReadList) {
        mappedReadList.forEach((read) -> methylPolarity += read.getMethylPolarity());
        this.mappedReadList.addAll(mappedReadList);
    }

    public String getId() {
        return id;
    }

    public List<Edge> getAdjEdges() {
        return adjEdges;
    }

    public List<MappedRead> getMappedReadList() {
        return mappedReadList;
    }

    public int getMethylPolarity() {
        return methylPolarity;
    }
}
