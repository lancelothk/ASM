package edu.cwru.cbc.ASM.detect.dataType;

import com.google.common.collect.Sets;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by kehu on 11/17/14.
 * Graph for detecting ASM.
 */
public class ASMGraph {
	private Set<Edge> edgeSet;
	private Map<String, Vertex> vertexMap; // id,read
	private Map<Integer, ClusterRefCpG> clusterRefCpGMap; // pos <--> count
	private Map<String, Vertex> clusterResult;// id,read
	// count how many tie situation occurs. For analysis use
	private int tieWeightCounter, tieIdCountCounter;
	private int originalEdgeCount;

	public ASMGraph(List<MappedRead> mappedReadList) {
		this.edgeSet = new LinkedHashSet<>();
		this.vertexMap = new HashMap<>();
		mappedReadList.forEach(r -> {
			if (vertexMap.put(r.getId(), new Vertex(r)) != null) {
				throw new RuntimeException("duplicate mapped read id!" + r.getId());
			}
		});
		// visit all possible edges once.
		for (int i = 0; i < mappedReadList.size(); i++) {
			for (int j = i + 1; j < mappedReadList.size(); j++) {
				double score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
				if (Double.compare(score, Double.MIN_VALUE) != 0) {
					edgeSet.add(new Edge(vertexMap.get(mappedReadList.get(i).getId()),
							vertexMap.get(mappedReadList.get(j).getId()), score));
				}
			}
		}

		this.originalEdgeCount = edgeSet.size();
	}

	private double checkCompatible(MappedRead readA, MappedRead readB) {
		if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
				readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
			int score = 0;
			for (CpG cpgA : readA.getCpgList()) {
				for (CpG cpgB : readB.getCpgList()) {
					if (cpgA.getPos() == cpgB.getPos()) {
						if (cpgA.getMethylStatus() != cpgB.getMethylStatus()) {
							score--;
						} else {
							score++;
						}
					}
				}
			}
			return score;
		} else {
			// don't have overlapped CpG, non-connected
			return Double.MIN_VALUE;
		}
	}

	public void cluster() {
		// merge vertexes connected by positive weight edge
		while (true) {
			// if edgeSet is empty
			if (edgeSet.size() == 0) {
				break;
			}

			List<Edge> maxEdgeList = getMaxWeightEdge(edgeSet);
			// if max weight < 0, stop merge.
			if (maxEdgeList.get(0).getWeight() < 0) {
				break;
			} else {
				// since edgeSet is not empty, getMaxWeightEdge always returns at least one element
				maxEdgeList = tieBreak(maxEdgeList);
				mergeVertex(maxEdgeList.get(0));
			}
		}

		setCoveredCpGMap();
		//refineClusterResult();
		//mergeNonOverlappedClusters();
		this.clusterResult = vertexMap;
	}

	private List<Edge> tieBreak(List<Edge> maxEdgeList) {
		if (maxEdgeList.size() != 1) {
			// multiple equal max weight edges.
			tieWeightCounter++;

			// First pick edge not connect to two clusters.
			// sort edge by id count. Smaller first
			maxEdgeList = getMinIdCountEdge(maxEdgeList);
			if (maxEdgeList.size() != 1) {
				tieIdCountCounter++;
			}
			//                    }
		}
		return maxEdgeList;
	}

	private List<Edge> getMaxWeightEdge(Set<Edge> itemList) {
		List<Edge> maxItemList = new ArrayList<>();
		double maxValue = Double.NEGATIVE_INFINITY;
		for (Edge item : itemList) {
			if (item.getWeight() > maxValue) {
				maxItemList.clear();
				maxItemList.add(item);
				maxValue = item.getWeight();
			} else if (item.getWeight() == maxValue) {
				maxItemList.add(item);
			}
		}
		if (maxItemList.size() == 0) {
			throw new RuntimeException("max item list is empty!");
		}
		return maxItemList;
	}

	private List<Edge> getMinIdCountEdge(List<Edge> itemList) {
		List<Edge> minItemList = new ArrayList<>();
		int minValue = itemList.get(0).getIdCount();
		for (Edge item : itemList) {
			if (item.getIdCount() < minValue) {
				minItemList.clear();
				minItemList.add(item);
				minValue = item.getIdCount();
			} else if (item.getIdCount() == minValue) {
				minItemList.add(item);
			}
		}
		return minItemList;
	}

	private void mergeVertex(Edge edge) {
		Vertex left = edge.getLeft();
		Vertex right = edge.getRight();
		// merge right to left vertex
		left.addMappedRead(right.getMappedReadList());
		left.addRefCpG(right.getRefCpGMap().values());
		//remove this edge and right vertex from graph
		edge.removeFromVertex();
		edgeSet.remove(edge);
		vertexMap.remove(right.getId());

		updateAdjEdges(left, right, edgeSet);

	}

	private void updateAdjEdges(Vertex left, Vertex right, Set<Edge> edgeSet) {
		HashMap<String, Edge> leftAdjEdgeMap = getStringEdgeHashMap(left);
		// update right vertex adj edges
		for (Edge rightEdge : right.getAdjEdgeList()) {
			// update new vertex
			rightEdge.replaceVertex(right, left);
			updateEdges(left, edgeSet, leftAdjEdgeMap, rightEdge);
		}
	}

	private void updateEdges(Vertex left, Set<Edge> edgeSet, HashMap<String, Edge> leftAdjEdgeMap, Edge rightEdge) {
		if (leftAdjEdgeMap.containsKey(rightEdge.getUniqueId())) {
			// update left edge weight
			Edge leftEdge = leftAdjEdgeMap.get(rightEdge.getUniqueId());
			leftEdge.setWeight((leftEdge.getWeight() + rightEdge.getWeight()) / 2);
			rightEdge.removeFromVertex();
			edgeSet.remove(rightEdge);
		} else {
			left.addEdge(rightEdge);
		}
	}

	private HashMap<String, Edge> getStringEdgeHashMap(Vertex left) {
		HashMap<String, Edge> leftAdjEdgeMap = new HashMap<>();
		left.getAdjEdgeList().forEach(e -> leftAdjEdgeMap.put(e.getUniqueId(), e));
		return leftAdjEdgeMap;
	}

	private void setCoveredCpGMap() {
		clusterRefCpGMap = new HashMap<>();
		for (Vertex vertex : vertexMap.values()) {
			for (RefCpG refCpG : vertex.getRefCpGMap().values()) {
				if (!clusterRefCpGMap.containsKey(refCpG.getPos())) {
					clusterRefCpGMap.put(refCpG.getPos(), new ClusterRefCpG(refCpG, vertex));
				} else {
					clusterRefCpGMap.get(refCpG.getPos()).addVertex(vertex);
				}
			}
		}
	}

	// merge until max cpg group coverage = 2
	private void refineClusterResult() {
		if (maxCpGGroupCoverage() > 2) {
			// merge vertexes connected by positive weight edge
			while (true) {
				// if edgeSet is empty
				if (edgeSet.size() == 0) {
					break;
				}

				List<Edge> maxEdgeList = getMaxWeightEdge(edgeSet);
				if (maxCpGGroupCoverage() == 2) {
					break;
				} else {
					// since edgeSet is not empty, getMaxWeightEdge always returns at least one element
					maxEdgeList = tieBreak(maxEdgeList);
					mergeVertex(maxEdgeList.get(0));
					setCoveredCpGMap();
				}
			}
		}
	}

	private void mergeNonOverlappedClusters() {
		while (true) {
			int prevSize = vertexMap.size();
			Set<Vertex> curr, next;
			List<ClusterRefCpG> clusterRefCpGList = clusterRefCpGMap.values()
					.stream()
					.filter(c -> c.getClusterCount() >= 2)
					.sorted()
					.collect(Collectors.toList());
			for (int i = 0; i < clusterRefCpGList.size() - 1; i++) {
				curr = clusterRefCpGList.get(i).getClusterSet();
				next = clusterRefCpGList.get(i + 1).getClusterSet();
				Set<Vertex> diff = Sets.symmetricDifference(curr, next);
				if (diff.size() == 2) {
					// merge two diff vertexes
					Vertex[] diffArray = diff.toArray(new Vertex[2]);
					mergeVertex(diffArray[0], diffArray[1]);
					break;
				} //else if (diff.size() > 2) {
				// more than 2 diff vertexes
				// impossible since after refineClusterResult(), no refCpG contains more than 2 vertexes.
				//}
				// else:
				// same set, do nothing or one set with 2 vertexes, one set with 1 vertex, do nothing
			}
			setCoveredCpGMap();
			// break loop when no more vertex can be merged.
			if (prevSize == vertexMap.size()) {
				break;
			}
		}
	}

	private int maxCpGGroupCoverage() {
		return clusterRefCpGMap.values()
				.stream()
				.max((c1, c2) -> c1.getClusterCount() - c2.getClusterCount())
				.get()
				.getClusterCount();
	}

	private void mergeVertex(Vertex a, Vertex b) {
		// merge b to a
		a.addMappedRead(b.getMappedReadList());
		a.addRefCpG(b.getRefCpGMap().values());
		vertexMap.remove(b.getId());

		// since only called in final merging, don't need to update edges.
	}

	public Map<Integer, ClusterRefCpG> getClusterRefCpGMap() {
		return clusterRefCpGMap;
	}

	public Map<String, Vertex> getClusterResult() {
		return clusterResult;
	}

	public int getTieWeightCounter() {
		return tieWeightCounter;
	}

	public int getTieIdCountCounter() {
		return tieIdCountCounter;
	}

	public int getOriginalEdgeCount() {
		return originalEdgeCount;
	}

	public double getNormMECSum() {
		return getMECSum() / getCpGSum();
	}

	public double getMECSum() {
		double mecSum = 0;
		for (Vertex vertex : clusterResult.values()) {
			mecSum += vertex.getMECScore();
		}
		return mecSum;
	}

	public int getCpGSum() {
		int cpgSum = 0;
		for (Vertex vertex : clusterResult.values()) {
			cpgSum += vertex.getCpGSum();
		}
		return cpgSum;
	}
}
