package edu.cwru.cbc.ASM.detect.dataType;

import edu.cwru.cbc.ASM.commons.methylation.RefCpG;

import javax.annotation.Nonnull;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by kehu on 12/8/14.
 * Used for record clusters in each RefCpG position.
 */
public class ClusterRefCpG implements Comparable<ClusterRefCpG> {
	private RefCpG refCpG;
	private Set<Vertex> clusterSet;


	public ClusterRefCpG(RefCpG refCpG, Vertex vertex) {
		this.refCpG = refCpG;
		this.clusterSet = new HashSet<>();
		this.clusterSet.add(vertex);
	}

	@Override
	public int compareTo(@Nonnull ClusterRefCpG o) {
		return this.refCpG.compareTo(o.refCpG);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (!(o instanceof ClusterRefCpG)) return false;

		ClusterRefCpG that = (ClusterRefCpG) o;
		return refCpG.equals(that.refCpG) && clusterSet.equals(that.clusterSet);

	}

	@Override
	public int hashCode() {
		int result = refCpG.hashCode();
		result = 31 * result + clusterSet.hashCode();
		return result;
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

	public RefCpG getRefCpG() {
		return refCpG;
	}
}
