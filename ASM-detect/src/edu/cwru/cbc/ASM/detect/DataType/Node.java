package edu.cwru.cbc.ASM.detect.DataType;

/**
 * Created by ke on 4/23/14.
 */
public class Node {
    private int index;
    private boolean isVisited;

    public Node(int index) {
        this.index = index;
        this.isVisited = false;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public boolean isVisited() {
        return isVisited;
    }

    public void setVisited(boolean isVisited) {
        this.isVisited = isVisited;
    }
}
