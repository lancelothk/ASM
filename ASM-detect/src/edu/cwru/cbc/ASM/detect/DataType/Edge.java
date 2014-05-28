package edu.cwru.cbc.ASM.detect.DataType;

/**
 * Created by lancelothk on 5/14/14.
 * Edge in graph.
 */
public class Edge {
    private Vertex left;
    private Vertex right;
    private int weight;

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

	public void setWeight(int weight) {
		this.weight = weight;
	}

	public boolean replaceVertex(Vertex o, Vertex n){
		if (o == n){
			throw new RuntimeException("old and new vertex are same!");
		}
		if (left == o){
			left = n;
			return true;
		}else if (right == o) {
			right = n;
			return true;
		}else {
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

	public int getIdCount(){
		return left.getIdList().size() + right.getIdList().size();
	}

	public void removeFromVertex() {
		this.left.removeEdge(this);
		this.right.removeEdge(this);
	}

	public String getUniqueId(){
		// smaller id first
		return left.getId() > right.getId()? String.format("%d-%d", left.getId(), right.getId())
				: String.format("%d-%d", right.getId(), left.getId());
	}
}
