package edu.cwru.cbc.ASM.detect;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.detect.DataType.CpGSite;
import edu.cwru.cbc.ASM.detect.DataType.Edge;
import edu.cwru.cbc.ASM.detect.DataType.MappedRead;
import edu.cwru.cbc.ASM.detect.DataType.Vertex;
import edu.cwru.cbc.ASM.detect.Utils.MappedReadFileLineProcessor;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Created by ke on 2/19/14.
 * ASM detection
 * Note:    1. mapped reads start pos is 1-based, end pos is 0-based.
 */
public class DetectASM implements Callable<String> {
    private List<Edge> edgeList;
	private Map<Long, Vertex> vertexMap;

    private File intervalFile;
    private long initPos;
    private BufferedWriter writer;
    private BufferedWriter group2Writer;

    public DetectASM(File intervalFile, long initPos, BufferedWriter writer, BufferedWriter group2Writer) {
        this.intervalFile = intervalFile;
        this.initPos = initPos;
        this.writer = writer;
        this.group2Writer = group2Writer;
    }

    public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
        long start = System.currentTimeMillis();
        // call on single file
//        detectASM.execute("/home/kehu/ASM_result/chr20-56897421-56898208.reads", 56897421);
//        detectASM.execute("ASM/testData/FindASM/test.reads", 1);
//		detectASM.execute(new File("/home/lancelothk/chr20_test/chr20-56895353-56895567"), 56895353);

        // call on folder
        // TODO mkdir if not exist
        String cellLine = "i90";
        String replicate = "r1";
        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(
                String.format("/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_summary", cellLine,
                              replicate)));
        summaryWriter.write("name\tlength\treadCount\tCpGCount\tGroupCount\n");
        BufferedWriter writer = new BufferedWriter(new FileWriter(
                String.format("/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_groups", cellLine,
                              replicate)));
        BufferedWriter group2Writer = new BufferedWriter(new FileWriter(
                String.format("/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_group2", cellLine,
                              replicate)));
        group2Writer.write("name\tlength\treadCount\tCpGCount\tGroupCount\t1stGroupSize\t2ndGroupSize\n");
        File path = new File(String.format("/home/kehu/experiments/ASM/result_%s_%s/intervals", cellLine, replicate));

//        ExecutorService executor = Executors.newCachedThreadPool();
        ExecutorService executor = Executors.newFixedThreadPool(4);
        List<Future<String>> futureList = new ArrayList<>();
        if (path.isDirectory()) {
			for (File file : path.listFiles()) {
				if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
						!file.getName().endsWith("intervalSummary")) {
//					if (file.getName().contains("22237970")){
						String[] items = file.getName().split("-");
                    Future<String> future = executor.submit(
                            new DetectASM(file, Long.parseLong(items[1]), writer, group2Writer));
                    futureList.add(future);
//					}
				}
			}
		}else {
			String[] items = path.getName().split("-");
            Future<String> future = executor.submit(
                    new DetectASM(path, Long.parseLong(items[1]), writer, group2Writer));
            futureList.add(future);
        }

        for (Future<String> stringFuture : futureList) {
            summaryWriter.write(stringFuture.get());
            System.out.println(stringFuture.get());
        }

		summaryWriter.close();
		writer.close();
		group2Writer.close();
		System.out.println(System.currentTimeMillis() - start + "ms");
	}


    @Override
    public String call() throws Exception {
        StringBuilder asm_result = new StringBuilder(intervalFile.getName() + "\n");
		StringBuilder summary = new StringBuilder(intervalFile.getName());
        String reference = readRef(intervalFile);
        List<CpGSite> cpgList = extractCpGSite(reference, initPos);
        List<MappedRead> mappedReadList = Files.asCharSource(intervalFile, Charsets.UTF_8).readLines(
                new MappedReadFileLineProcessor());
        associateReadWithCpG(cpgList, mappedReadList);
		String[] items = intervalFile.getName().split("-");
		if (items.length != 3){
            throw new RuntimeException("illegal interval file name format!");
        }
		summary.append("\t" + (Integer.parseInt(items[2]) - Integer.parseInt(items[1]) + 1));
		summary.append("\t" + mappedReadList.size());
		summary.append("\t" + cpgList.size());

		vertexMap = new HashMap<>();
		edgeList = new ArrayList<>();
        constructGraph(vertexMap, edgeList, mappedReadList);
		getClusters(asm_result);

		int groupCount = 0;

//        int minGroupSize = Integer.MIN_VALUE;
		for (Vertex vertex : vertexMap.values()) {
//            minGroupSize = minGroupSize > vertex.getIdList().size()? vertex.getIdList().size(): minGroupSize;
			groupCount++;
			for (Long id : vertex.getIdList()) {
                asm_result.append(id + ",");
			}
            asm_result.append("\n");
		}
        asm_result.append("Number of groups:\t" + groupCount + "\n");
		summary.append("\t" + groupCount + "\n");
		String intervalSummary = summary.toString();
//        if (minGroupSize < 4){
//            return "";
//        }else {
            if (vertexMap.values().size() == 2){
                group2Writer.write(intervalSummary + "\t");
                for (Vertex vertex : vertexMap.values()) {
                    group2Writer.write(vertex.getIdList().size() + "\t");
                }
                group2Writer.write("\n");
            }
            writer.write(asm_result.toString());
            return summary.toString();
//        }
    }

    private void getClusters(StringBuilder asm_result) throws IOException {
        long start = System.currentTimeMillis();
        // count how many tie situation occurs. For analysis use
		int tieWeightCounter = 0, tieIdCountCounter = 0;
		// merge vertexes connected by positive weight edge
		while (true) {
            if ((System.currentTimeMillis() - start) > 30 * 1000) {
                System.out.println(intervalFile.getName());
            }
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
				mergeVertex(edge);
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
				mergeVertex(maxEdgeList.get(0));
				// record the frequency of tied id count. Maybe solve by considering inner cluster weight
				// or some other properties
				if (maxEdgeList.get(0).getIdCount() == maxEdgeList.get(1).getIdCount()){
					tieIdCountCounter++;
				}
			}
		}
        asm_result.append(String.format("tie weight counter:%d\n", tieWeightCounter));
        asm_result.append(String.format("tie id counter:%d\n", tieIdCountCounter));
	}

	private void mergeVertex(Edge edge){
		Vertex left = edge.getLeft();
		Vertex right = edge.getRight();
		// merge right to left vertex
		left.addIds(right.getIdList());
		left.addInnerWeight(edge.getWeight());
		//remove this edge and right vertex from graph
		edge.removeFromVertex();
		edgeList.remove(edge);
		vertexMap.remove(right.getId());
		// update right vertex adj edges
		for (Edge edgeR : right.getAdjEdges()) {
			// update new vertex
			if (!edgeR.replaceVertex(right, left)) {
				throw new RuntimeException("fail to update new vertex in adj edge!");
			}
		}
		// update left vertex with new edges
		for (Edge adjEdge : right.getAdjEdges()) {
			left.addEdge(adjEdge);
		}
		// update edges of vertex which connect both left and right
		updateAndRemoveDupEdge();
	}

	private void updateAndRemoveDupEdge(){
		Map<String, Edge> edgeMap = new HashMap<>();
		Iterator<Edge> edgeIterator = edgeList.iterator();
		while (edgeIterator.hasNext()){
			Edge edgeB = edgeIterator.next();
			// duplicate edge
			if (edgeMap.containsKey(edgeB.getUniqueId())){
				Edge edgeA = edgeMap.get(edgeB.getUniqueId());
				// keep edgeA, remove edgeB.
				edgeA.setWeight(edgeA.getWeight() + edgeB.getWeight());
				edgeB.removeFromVertex();
				edgeIterator.remove(); // remove from edgelist
			}else {
				edgeMap.put(edgeB.getUniqueId(), edgeB);
			}
		}
		edgeList = new ArrayList<>(edgeMap.values());
	}

	/**
	 * select edge/edges with max weight in the list.
     * @param edgeList list contains edges
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

    private String readRef(File intervalFile) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(intervalFile));
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
