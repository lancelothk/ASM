package edu.cwru.cbc.ASM.detect.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lancelothk on 5/14/14.
 * Vertex in graph.
 */
public class Vertex {
    private long id;
    private List<Edge> adjEdges;
    private List<Long> idList;
	private int innerWeightSum = 0;

    public Vertex(long id) {
        this.id = id;
        this.adjEdges = new ArrayList<>();
        this.idList = new ArrayList<>();
		this.idList.add(id);
    }

    public void addEdge(Edge edge){
        this.adjEdges.add(edge);
    }

	public void removeEdge(Edge edge){
		this.adjEdges.remove(edge);
	}

	public void addIds(List<Long> idList){
		this.idList.addAll(idList);
	}

	public long getId() {
		return id;
	}

	public void addInnerWeight(int weight){
		this.innerWeightSum += weight;
	}

	public List<Edge> getAdjEdges() {
		return adjEdges;
	}

	public List<Long> getIdList() {
		return idList;
	}
}
