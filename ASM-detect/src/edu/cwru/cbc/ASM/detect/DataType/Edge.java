package edu.cwru.cbc.ASM.detect.DataType;

/**
 * Created by lancelothk on 5/14/14.
 * Edge in graph.
 */
public class Edge {
    private Vertex left;
    private Vertex right;
    private int weight;
    private String id = null;

    public Edge(Vertex left, Vertex right, int weight) {
		if (left == right){
			throw new RuntimeException("same vertex in one edge!");
		}
        left.addEdge(this);
        right.addEdge(this);
        this.left = left;
        this.right = right;
        this.weight = weight;
    }

    public boolean replaceVertex(Vertex o, Vertex n) {
        if (o == n) {
            throw new RuntimeException("old and new vertex are same!");
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

    public int getWeight() {
        return weight;
    }

    public void setWeight(int weight) {
        this.weight = weight;
    }

	public int getIdCount(){
//		return Math.min(left.getIdList().size(), right.getIdList().size());
		return left.getIdList().size() + right.getIdList().size();
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
        this.id =
                left.getId() <= right.getId() ? left.getId() + "-" + right.getId() : right.getId() + "-" + left.getId();
        return id;
    }
}
