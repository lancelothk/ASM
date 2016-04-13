package edu.cwru.cbc.ASM.detect.dataType;

/**
 * Created by lancelothk on 5/14/14.
 * Edge in graph.
 */
public class Edge {
	private Vertex left;
	private Vertex right;
	private double weight;
	private String id;

	public Edge(Vertex left, Vertex right, double weight) {
		if (left == right) {
			throw new RuntimeException("same vertex in one edge!" + left.getId());
		}
		this.left = left;
		this.right = right;
		this.weight = weight;
		left.addEdge(this);
		right.addEdge(this);
		calcId();
	}

	private String calcId() {
		// smaller id first
		this.id = left.getId().compareTo(right.getId()) < 0 ?
				left.getId() + "-" + right.getId() : right.getId() + "-" + left.getId();
		return id;
	}

	public void replaceVertex(Vertex oldVertex, Vertex newVertex) {
		if (oldVertex == newVertex) {
			throw new RuntimeException("oldVertex and new vertex are same!" + oldVertex.getId());
		}
		if (left == oldVertex) {
			left = newVertex;
			calcId();
		} else if (right == oldVertex) {
			right = newVertex;
			calcId();
		} else {
			throw new RuntimeException("This edge do not contain oldVertex!\t" + oldVertex.getId());
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

	public int getIdCount() {
		return left.getMappedReadList().size() + right.getMappedReadList().size();
	}

	public void removeFromVertex() {
		this.left.removeEdge(this);
		this.right.removeEdge(this);
	}

	public String getUniqueId() {
		return id == null ? calcId() : id;
	}
}