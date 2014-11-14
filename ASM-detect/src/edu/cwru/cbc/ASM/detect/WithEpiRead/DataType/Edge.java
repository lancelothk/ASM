package edu.cwru.cbc.ASM.detect.WithEpiRead.DataType;

/**
 * Created by kehu on 11/13/14.
 * Edge for Matrxi represented graph
 */
public class Edge {
    private int left;
    private int right;
    private int weight;

    public Edge(int left, int right, int weight) {
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

    public int getWeight() {
        return weight;
    }
}
