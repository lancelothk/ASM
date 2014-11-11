package edu.cwru.cbc.ASM.detect.DataType;

import java.util.*;

/**
 * Created by lancelothk on 5/14/14.
 * Vertex in graph.
 */
public class Vertex {
    private long id;
    private List<Edge> adjEdges;
    private List<MappedRead> mappedReadList;
    private Map<Long, CpGSite> cpGSiteMap;
    private int methylPolarity;
    private int innerWeightSum;

    public Vertex(MappedRead mappedRead) {
        this.id = mappedRead.getId();
        this.adjEdges = new ArrayList<>();
        this.mappedReadList = new ArrayList<>();
        this.mappedReadList.add(mappedRead);
        this.cpGSiteMap = new HashMap<>();
        this.methylPolarity = mappedRead.getMethylPolarity();
        addCpGSites(mappedRead.getCpgList());
    }

    public void addCpGSites(Collection<CpGSite> cpGSiteList) {
        for (CpGSite cpGSite : cpGSiteList) {
            if (!cpGSiteMap.containsKey(cpGSite.getPos())) {
                cpGSiteMap.put(cpGSite.getPos(), cpGSite);
            }
        }
    }

    public Collection<CpGSite> getCoveredCpGSites() {
        return cpGSiteMap.values();
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

    public long getId() {
        return id;
    }

    public void addInnerWeight(int weight) {
        this.innerWeightSum += weight;
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
