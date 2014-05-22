package ASM.DataType;

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

    public Vertex(long id) {
        this.id = id;
        this.adjEdges = new ArrayList<>();
        this.idList = new ArrayList<>();
    }

    public void addEdge(Edge edge){
        this.adjEdges.add(edge);
    }
}
