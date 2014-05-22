package ASM.DataType;

/**
 * Created by lancelothk on 5/14/14.
 * Edge in graph.
 */
public class Edge {
    private Vertex left;
    private Vertex right;
    private int weight;

    public Edge(Vertex left, Vertex right, int weight) {
        left.addEdge(this);
        right.addEdge(this);
        this.left = left;
        this.right = right;
        this.weight = weight;
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
}
