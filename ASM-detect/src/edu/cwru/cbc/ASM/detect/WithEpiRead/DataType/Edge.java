package edu.cwru.cbc.ASM.detect.WithEpiRead.DataType;

/**
 * Created by kehu on 11/13/14.
 * Edge for Matrxi represented graph
 */
public class Edge implements Comparable<Edge> {
    private int left;
    private int right;
    private double weight;

    public Edge(int left, int right, double weight) {
        this.left = left;
        this.right = right;
        this.weight = weight;
    }

    public int getLeft() {
        return left;
    }

    public int getRight() {
        return right;
    }

    public double getWeight() {
        return weight;
    }

    @Override
    public int compareTo(Edge o) {
        // larger first
        return Double.compare(o.getWeight(), this.weight);
    }
}
