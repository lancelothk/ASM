package ASM.Execution;

import ASM.DataType.*;
import ASM.Utils.MappedReadFileLineProcessor;
import com.google.common.base.Charsets;
import com.google.common.io.Files;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by ke on 2/19/14.
 * ASM detection
 * Note:    1. mapped reads start pos is 1-based, end pos is 0-based.
 */
public class FindASM {
    private enum Compatibility {
        COMPATIBLE, NONCOMPATIBLE, NONCONNECTED, SELF
    }

    private Compatibility[][] matrix;

    public static void main(String[] args) throws IOException {
        FindASM findASM = new FindASM();
//        findASM.execute("/home/kehu/ASM_result/chr20-56897421-56898208.reads", 56897421);
//        findASM.execute("ASM/testData/FindASM/test.reads", 1);
		findASM.execute("/home/lancelothk/chr20-56897421-56898208_short.reads", 56897421);
    }

    public void execute(String intervalFileName, long initPos) throws IOException {
        String reference = readRef(intervalFileName);
        List<CpGSite> cpgList = extractCpGSite(reference, initPos);
        List<MappedRead> mappedReadList = Files.asCharSource(new File(intervalFileName), Charsets.UTF_8).readLines(
                new MappedReadFileLineProcessor());
        associateReadWithCpG(cpgList, mappedReadList);
        for (int i = 0; i < mappedReadList.size(); i++) {
            if (mappedReadList.get(i).getCpgList() == null) {
                mappedReadList.remove(i--);
            }
        }
        //Collections.sort(mappedReadList, new MappedReadComparaterByCpG());
		Map<Long, Vertex> vertexMap = new HashMap<>();
		List<Edge> edgeList = new ArrayList<>();
        constructGraph(vertexMap, edgeList, mappedReadList);
		getClusters(vertexMap, edgeList);
		for (Vertex vertex : vertexMap.values()) {
			for (Long id : vertex.getIdList()) {
				System.out.print(id + ",");
			}
			System.out.println("");
		}

//        List<List<Node>> resultList = detectASM(mappedReadList);
//        for (List<Node> nodes : resultList) {
//            System.out.printf("size:\t%d\t->", nodes.size());
//            for (Node node : nodes) {
//                System.out.printf("%d\t", mappedReadList.get(node.getIndex()).getId());
//            }
//            System.out.println();
//        }
    }

    private void getClusters(Map<Long, Vertex> vertexMap, List<Edge> edgeList){
		// count how many tie situation occurs. For analysis use
		int tieWeightCounter = 0, tieIdCountCounter = 0;
		// merge vertexes connected by positive weight edge
		while (true) {
			// if edgeList is empty
			if (edgeList.size() == 0){
				break;
			}
			List<Edge> maxEdgeList = getMaxEdge(edgeList);
			// if max weight <= 0, stop merge.
			if (maxEdgeList.get(0).getWeight() <= 0) {
				break;
			} else if (maxEdgeList.size() == 1) {
				// unique max weight, merge two vertex on this edge and updating adj edges.
				Edge edge = maxEdgeList.get(0);
				mergeVertex(edge, vertexMap, edgeList);
			} else {
				// multiple equal max weight edges. First pick edge not connect to two clusters.
				tieWeightCounter++;
				// sort edge by id count
				Collections.sort(maxEdgeList, new Comparator<Edge>() {
					@Override
					public int compare(Edge o1, Edge o2) {
						return o1.getIdCount() - o2.getIdCount();
					}
				});
				// pick min id count edge as merge candidate. i.e. prefer smaller cluster.
				mergeVertex(maxEdgeList.get(0), vertexMap ,edgeList);
				// record the frequency of tied id count. Maybe solve by considering inner cluster weight
				// or some other properties
				if (maxEdgeList.get(0).getIdCount() == maxEdgeList.get(1).getIdCount()){
					tieIdCountCounter++;
				}
			}
		}
		System.out.printf("tie weight counter:%d\n", tieWeightCounter);
		System.out.printf("tie id counter:%d\n", tieIdCountCounter);
	}

	private void mergeVertex(Edge edge, Map<Long, Vertex> vertexMap, List<Edge> edgeList){
		Vertex left = edge.getLeft();
		Vertex right = edge.getRight();
		// merge right to left vertex
		left.addIds(right.getIdList());
		//remove this edge and right vertex from graph
		edge.removeFromVertex();
		edgeList.remove(edge);
		vertexMap.remove(right.getId());
		// update right vertex adj edges
		for (Edge edgeR : right.getAdjEdges()) {
			// update new vertex
			if (!edgeR.replaceVertex(right, left)) {
//				throw new RuntimeException("fail to update new vertex in adj edge!");
			}else {
				for (Edge adjEdge : right.getAdjEdges()) {
					left.addEdge(adjEdge);
				}
			}
			// update new weight
			edgeR.setWeight(edgeR.getWeight() + edge.getWeight());
		}
		// update edges of vertex which connect both left and right
		updateAndRemoveDupEdge(edgeList, edge);
	}

	private void updateAndRemoveDupEdge(List<Edge> edgeList, Edge targetEdge){
		Map<String, Edge> edgeMap = new HashMap<>();
		Iterator<Edge> edgeIterator = edgeList.iterator();
		while (edgeIterator.hasNext()){
			Edge edge = edgeIterator.next();
			// duplicate edge
			if (edgeMap.containsKey(edge.getUniqueId())){
				Edge prevEdge = edgeMap.get(edge.getUniqueId());
				// keep edgeA, remove edgeB. edge.weight should be excluded here.
				prevEdge.setWeight(prevEdge.getWeight() + edge.getWeight() - targetEdge.getWeight());
				edge.removeFromVertex();
				edgeIterator.remove();
			}else {
				edgeMap.put(edge.getUniqueId(), edge);
			}
		}
		edgeList = new ArrayList<>(edgeMap.values());
	}

	/**
	 * select edge/edges with max weight in the list.
	 * @param edgeList
	 * @return	list contains edge/edges with max weight
	 */
	private List<Edge> getMaxEdge(List<Edge> edgeList){
		List<Edge> maxEdgeList = new ArrayList<>();
		int maxWeight = Integer.MIN_VALUE;
		for (Edge edge : edgeList) {
			if (edge.getWeight() > maxWeight){
				maxEdgeList.clear();
				maxEdgeList.add(edge);
				maxWeight = edge.getWeight();
			}else if (edge.getWeight() == maxWeight){
				maxEdgeList.add(edge);
			}
		}
		return maxEdgeList;
	}

    private List<List<Node>> detectASM(List<MappedRead> mappedReadList) {
        List<List<Node>> resultList = new ArrayList<>();
        List<Node> nodeList = new ArrayList<>();
        for (int i = 0; i < mappedReadList.size(); i++) {
            nodeList.add(new Node(i));
        }
        while (true) {
            List<Node> result = getGroup(nodeList);
            if (result == null) {
                break;
            } else {
                resultList.add(result);
                for (Node node : result) {
                    node.setVisited(true);
                }
            }
        }
        return resultList;
    }

    private List<Node> getGroup(List<Node> nodeList) {
        List<Node> result = new ArrayList<>();
        List<Node> seedList = new ArrayList<>();
        Node firstNode = null;
        for (Node aNodeList : nodeList) {
            if (!aNodeList.isVisited()) {
                firstNode = aNodeList;
                break;
            }
        }
        if (firstNode != null) {
            seedList.add(firstNode);
            result.add(firstNode);
        } else {
            // all node are visited
            return null;
        }
        while (true) {
            if (seedList.size() == 0) {
                break;
            } else {
                addCandidates(nodeList, result, seedList, seedList.get(0));
            }
        }
        return result;
    }

    private void addCandidates(List<Node> nodeList, List<Node> result, List<Node> seedList, Node seed) {
        // check all other reads
        for (int i = 0; i < matrix[seed.getIndex()].length; i++) {
            if (!nodeList.get(i).isVisited() && !result.contains(nodeList.get(i)) &&
                    !seedList.contains(nodeList.get(i)) && matrix[seed.getIndex()][i] == Compatibility.COMPATIBLE) {
                //if candidate has no contradiction with result member, add to result
                if (!hasContradiction(result, i)) {
                    seedList.add(nodeList.get(i));
                    result.add(nodeList.get(i));
                }
            }
        }
        seedList.remove(seed);
    }

    private boolean hasContradiction(List<Node> result, int candidate) {
        for (Node member : result) {
            if (matrix[candidate][member.getIndex()] == Compatibility.NONCOMPATIBLE) {
                return true;
            }
        }
        return false;
    }

    private void constructGraph(Map<Long, Vertex> vertexMap, List<Edge> edgeList, List<MappedRead> mappedReadList) {
        // visit all possible edges once.
        for (int i = 0; i < mappedReadList.size(); i++) {
            for (int j = i+1; j < mappedReadList.size(); j++) {
                int score = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
                if (score != Integer.MIN_VALUE){
                    long idI = mappedReadList.get(i).getId(), idJ = mappedReadList.get(j).getId();
                    if (!vertexMap.containsKey(idI)){
                        vertexMap.put(idI, new Vertex(idI));
                    }
                    if (!vertexMap.containsKey(idJ)){
                        vertexMap.put(idJ, new Vertex(idJ));
                    }
                    edgeList.add(new Edge(vertexMap.get(idI), vertexMap.get(idJ), score));
                }
            }
        }
    }

    private int checkCompatible(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            int score = 0;
            for (CpGSite cpgA : readA.getCpgList()) {
                for (CpGSite cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos()) {
                        if (cpgA.isMethylated() != cpgB.isMethylated()) {
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
            return Integer.MIN_VALUE;
        }
    }

    private String readRef(String intervalFileName) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(intervalFileName));
        String line = bufferedReader.readLine();
        String[] items = line.split("\t");
        if (items.length != 2) {
            throw new RuntimeException("invalid reference line in interval read file!");
        } else {
            return items[1].toUpperCase();
        }
    }

    private void associateReadWithCpG(List<CpGSite> cpgList, List<MappedRead> mappedReadList) {
        for (MappedRead read : mappedReadList) {
            for (CpGSite cpg : cpgList) {
                if (cpg.getPos() >= read.getStart() && cpg.getPos() + 1 <= read.getEnd()) {
                    cpg.addAssociatedRead(read);
                    read.addCpG(new CpGSite(cpg.getPos(), read.getCpGMethylStatus(cpg.getPos())));
                }
            }
        }
    }

    private List<CpGSite> extractCpGSite(String reference, long initPos) {
        reference = reference.replace(" ", "");
        List<CpGSite> cpgList = new ArrayList<>();
        for (int i = 0; i < reference.length() - 1; i++) {
            if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
                cpgList.add(new CpGSite(initPos + i));
            }
        }
        return cpgList;
    }
}
