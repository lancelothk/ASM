package edu.cwru.cbc.ASM.detect.dataType;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by kehu on 12/8/14.
 * Used for record clusters in each RefCpG position.
 */
public class ClusterRefCpG implements Comparable<ClusterRefCpG> {
	private int pos;
	private Set<Vertex> clusterSet;


	public ClusterRefCpG(int pos, Vertex vertex) {
		this.pos = pos;
		this.clusterSet = new HashSet<>();
		this.clusterSet.add(vertex);
	}

	@Override
	public int compareTo(ClusterRefCpG o) {
		return this.pos - o.pos;
	}

	public int getClusterCount() {
		return clusterSet.size();
	}

	public Set<Vertex> getClusterSet() {
		return clusterSet;
	}

	public void addVertex(Vertex vertex) {
		clusterSet.add(vertex);
	}
}
