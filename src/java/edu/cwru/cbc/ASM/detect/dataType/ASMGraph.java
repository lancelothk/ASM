package edu.cwru.cbc.ASM.detect.dataType;

import com.google.common.collect.Sets;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by kehu on 11/17/14.
 * Graph for detecting ASM.
 */
public class ASMGraph {
	private List<Edge> edgeList;
	private Map<Integer, Vertex> vertexMap;
	private Map<Integer, ClusterRefCpG> clusterRefCpGMap; // pos <--> count
	private Map<Integer, Vertex> clusterResult;
	// count how many tie situation occurs. For analysis use
	private int tieWeightCounter, tieIdCountCounter;
	private int originalVertexCount, originalEdgeCount;

	public ASMGraph(List<MappedRead> mappedReadList) {
		this.edgeList = new ArrayList<>();
		this.vertexMap = new HashMap<>();
		// visit all possible edges once.
		for (int i = 0; i < mappedReadList.size(); i++) {
			for (int j = i + 1; j < mappedReadList.size(); j++) {
				double score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
				if (Double.compare(score, Double.MIN_VALUE) != 0) {
					int idI = mappedReadList.get(i).getId(), idJ = mappedReadList.get(j).getId();
					if (!vertexMap.containsKey(idI)) {
						vertexMap.put(idI, new Vertex(mappedReadList.get(i)));
					}
					if (!vertexMap.containsKey(idJ)) {

						vertexMap.put(idJ, new Vertex(mappedReadList.get(j)));
					}
					edgeList.add(new Edge(vertexMap.get(idI), vertexMap.get(idJ), score));
				}
			}
		}
		this.originalVertexCount = vertexMap.size();
		this.originalEdgeCount = edgeList.size();
	}

	private double checkCompatible(MappedRead readA, MappedRead readB) {
		try {
			if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
					readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
				int score = 0;
				for (CpG cpgA : readA.getCpgList()) {
					for (CpG cpgB : readB.getCpgList()) {
						if (cpgA.getPos() == cpgB.getPos()) {
							if (cpgA.getMethylStatus() != cpgB.getMethylStatus()) {
								assert cpgA.getMethylStatus() != MethylStatus.N;
								assert cpgB.getMethylStatus() != MethylStatus.N;
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
		} catch (Exception e) {
			System.out.printf("");
			throw e;
		}
	}

	public void cluster() throws IOException {
		// merge vertexes connected by positive weight edge
		while (true) {
			// if edgeList is empty
			if (edgeList.size() == 0) {
				break;
			}

			List<Edge> maxEdgeList = getMaxFromList(edgeList, Edge::getWeight);
			// if max weight <= 0, stop merge.
			if (maxEdgeList.get(0).getWeight() < 0) {
				break;
			} else {
				// since edgeList is not empty, getMaxFromList always returns at least one element
				if (maxEdgeList.size() != 1) {
					// multiple equal max weight edges.
					tieWeightCounter++;

					// First pick edge not connect to two clusters.
					// sort edge by id count. Smaller first
					maxEdgeList = getMinFromList(maxEdgeList, Edge::getIdCount);
					if (maxEdgeList.size() != 1) {
						tieIdCountCounter++;
					}
					//                    }
				}
				mergeVertex(maxEdgeList.get(0));
			}
		}

		setCoveredCpGMap();
		refineClusterResult();
		mergeNonoverlappedClusters();
		this.clusterResult = vertexMap;
	}

	/**
	 * select items with max property value in the list.
	 *
	 * @return list contains items with max weight
	 */
	private <T> List<T> getMaxFromList(List<T> itemList, Function<T, Double> getProperty) {
		List<T> maxItemList = new ArrayList<>();
		double maxValue = getProperty.apply(itemList.get(0));
		for (T item : itemList) {
			if (getProperty.apply(item) > maxValue) {
				maxItemList.clear();
				maxItemList.add(item);
				maxValue = getProperty.apply(item);
			} else if (getProperty.apply(item) == maxValue) {
				maxItemList.add(item);
			}
		}
		if (maxItemList.size() == 0) {
			throw new RuntimeException("max item list is empty!");
		}
		return maxItemList;
	}

	/**
	 * select items with min property value in the list.
	 *
	 * @return list contains items with min weight
	 */
	private <T> List<T> getMinFromList(List<T> itemList, Function<T, Integer> getProperty) {
		List<T> minItemList = new ArrayList<>();
		int minValue = getProperty.apply(itemList.get(0));
		for (T item : itemList) {
			if (getProperty.apply(item) < minValue) {
				minItemList.clear();
				minItemList.add(item);
				minValue = getProperty.apply(item);
			} else if (getProperty.apply(item) == minValue) {
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
		edgeList.remove(edge);
		vertexMap.remove(right.getId());
		// update right vertex adj edges
		for (Edge edgeR : right.getAdjEdges()) {
			// update new vertex
			edgeR.replaceVertex(right, left);
		}
		// update left vertex with new edges
		right.getAdjEdges().forEach(left::addEdge);
		// update edges of vertex which connect both left and right
		updateAndRemoveDupEdge(edgeList);
	}

	private void setCoveredCpGMap() {
		clusterRefCpGMap = new HashMap<>();
		for (Vertex vertex : vertexMap.values()) {
			for (RefCpG refCpG : vertex.getRefCpGMap().values()) {
				if (!clusterRefCpGMap.containsKey(refCpG.getPos())) {
					clusterRefCpGMap.put(refCpG.getPos(), new ClusterRefCpG(refCpG.getPos(), vertex));
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
				// if edgeList is empty
				if (edgeList.size() == 0) {
					break;
				}

				List<Edge> maxEdgeList = getMaxFromList(edgeList, Edge::getWeight);
				// if max weight <= 0, stop merge.
				if (maxCpGGroupCoverage() == 2) {
					break;
				} else {
					// since edgeList is not empty, getMaxFromList always returns at least one element
					if (maxEdgeList.size() != 1) {
						// multiple equal max weight edges.
						tieWeightCounter++;

						// First pick edge not connect to two clusters.
						// sort edge by id count. Smaller first
						maxEdgeList = getMinFromList(maxEdgeList, Edge::getIdCount);
						if (maxEdgeList.size() != 1) {
							tieIdCountCounter++;
						}
						//                    }
					}
					mergeVertex(maxEdgeList.get(0));
					setCoveredCpGMap();
				}
			}
		}
	}

	private void mergeNonoverlappedClusters() {
		while (true) {
			int prevSize = vertexMap.size();
			Set<Vertex> curr, next;
			List<ClusterRefCpG> clusterRefCpGList = clusterRefCpGMap.values()
					.stream()
					.filter(c -> c.getClusterCount() >= 2)
					.collect(Collectors.toList());
			clusterRefCpGList.sort(ClusterRefCpG::compareTo);
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
				// TODO how to merge?
				//}
				// else:
				// same set, do nothing
				// or
				// one set with 2 vertexes, one set with 1 vertex, do nothing
			}
			setCoveredCpGMap();
			// break loop when no more vertex can be merged.
			if (prevSize == vertexMap.size()) {
				break;
			}
		}
	}

	private void updateAndRemoveDupEdge(List<Edge> edgeList) {
		Map<String, Edge> edgeMap = new HashMap<>();
		Iterator<Edge> edgeIterator = edgeList.iterator();
		while (edgeIterator.hasNext()) {
			Edge edgeB = edgeIterator.next();
			// duplicate edge
			if (edgeMap.containsKey(edgeB.getUniqueId())) {
				Edge edgeA = edgeMap.get(edgeB.getUniqueId());
				// keep edgeA, remove edgeB.
				edgeA.setWeight(edgeA.getWeight() + edgeB.getWeight());
				edgeB.removeFromVertex();
				edgeIterator.remove(); // remove from edgelist
			} else {
				edgeMap.put(edgeB.getUniqueId(), edgeB);
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

	public Map<Integer, Vertex> getClusterResult() {
		return clusterResult;
	}

	public int getTieWeightCounter() {
		return tieWeightCounter;
	}

	public int getTieIdCountCounter() {
		return tieIdCountCounter;
	}

	public int getOriginalVertexCount() {
		return originalVertexCount;
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
