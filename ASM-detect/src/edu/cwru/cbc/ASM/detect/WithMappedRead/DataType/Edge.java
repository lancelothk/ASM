package edu.cwru.cbc.ASM.detect.WithMappedRead.DataType;

import static java.lang.Math.abs;

/**
 * Created by lancelothk on 5/14/14.
 * Edge in graph.
 */
public class Edge {
    private Vertex left;
    private Vertex right;
    private double weight;
    private String id = null;

    public Edge(Vertex left, Vertex right, double weight) {
        if (left == right){
            throw new RuntimeException("same vertex in one edge!" + left.getId());
        }
        left.addEdge(this);
        right.addEdge(this);
        this.left = left;
        this.right = right;
        this.weight = weight;
        calcId();
    }

    public boolean replaceVertex(Vertex o, Vertex n) {
        if (o == n) {
            throw new RuntimeException("old and new vertex are same!" + o.getId());
        }
        if (left == o) {
            left = n;
            calcId();
            return true;
        } else if (right == o) {
            right = n;
            calcId();
            return true;
        } else {
            return false;
        }
    }

    public Vertex getLeft() {
        return left;
    }

    public Vertex getRight() {
        return right;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

	public int getIdCount(){
        return left.getMappedReadList().size() + right.getMappedReadList().size();
    }

    public int getMethylPolarityAbs() {
        return abs(left.getMethylPolarity() + right.getMethylPolarity());
    }

	public void removeFromVertex() {
		this.left.removeEdge(this);
		this.right.removeEdge(this);
	}

	public String getUniqueId(){
        return id == null ? calcId() : id;
    }

    private String calcId() {
        // smaller id first
        this.id = left.getId().compareTo(right.getId()) < 0 ?
                left.getId() + "-" + right.getId() : right.getId() + "-" + left.getId();
        return id;
    }
}