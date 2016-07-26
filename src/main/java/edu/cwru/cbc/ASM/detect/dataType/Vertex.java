package edu.cwru.cbc.ASM.detect.dataType;

import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.util.*;

/**
 * Created by lancelothk on 5/14/14.
 * Vertex in graph.
 */
public class Vertex {
	private String id;
	private List<Edge> adjEdgeList;
	private List<MappedRead> mappedReadList;
	private Map<Integer, RefCpG> refCpGMap; // position-RefCpG

	public Vertex(MappedRead mappedRead) {
		this.id = mappedRead.getId();
		this.adjEdgeList = new ArrayList<>();
		this.mappedReadList = new ArrayList<>();
		this.mappedReadList.add(mappedRead);
		this.refCpGMap = new HashMap<>();
		this.addCpG(mappedRead.getCpgList());
	}

	public void addMappedRead(MappedRead mappedRead) {
		if (mappedRead.getId().equals(this.id)) {
			throw new RuntimeException("duplicate mapped read id!" + this.id);
		}
		this.mappedReadList.add(mappedRead);
		this.addCpG(mappedRead.getCpgList());
	}

	public void addCpG(Collection<CpG> cpgList) {
		cpgList.stream()
				.filter(cpg -> !refCpGMap.containsKey(cpg.getPos()))
				.forEach(cpg -> {
					RefCpG refCpG = new RefCpG(cpg.getPos());
					refCpG.addCpG(cpg);
					refCpGMap.put(cpg.getPos(), refCpG);
				});
	}

	public void addRefCpG(Collection<RefCpG> refCpgList) {
		for (RefCpG refCpG : refCpgList) {
			if (refCpGMap.containsKey(refCpG.getPos())) {
				refCpGMap.get(refCpG.getPos()).addMethylCount(refCpG.getMethylCount());
				refCpGMap.get(refCpG.getPos()).addNonMethylCount(refCpG.getCoveredCount() - refCpG.getMethylCount());
			} else {
				refCpGMap.put(refCpG.getPos(), refCpG);
			}
		}
	}

	public Map<Integer, RefCpG> getRefCpGMap() {
		return refCpGMap;
	}

	public void addEdge(Edge edge) {
		this.adjEdgeList.add(edge);
	}

	public void removeEdge(Edge edge) {
		this.adjEdgeList.remove(edge);
	}

	public void addMappedReads(List<MappedRead> mappedReadList) {
		for (MappedRead mappedRead : mappedReadList) {
			this.addMappedRead(mappedRead);
		}
	}

	public String getId() {
		return id;
	}

	public List<Edge> getAdjEdgeList() {
		return adjEdgeList;
	}

	public List<MappedRead> getMappedReadList() {
		return mappedReadList;
	}

	public double getMECScore() {
		double mec = 0;
		for (MappedRead mappedRead : mappedReadList) {
			for (CpG cpG : mappedRead.getCpgList()) {
				MethylStatus refMethylStatus = refCpGMap.get(cpG.getPos()).getMajorMethylStatus();
				assert cpG.getMethylStatus() != MethylStatus.N;
				assert refMethylStatus != MethylStatus.N;
				if (refMethylStatus == MethylStatus.E) {
					mec += 0.5;
				} else if (cpG.getMethylStatus() != refMethylStatus && cpG.getMethylStatus() != MethylStatus.E) {
					mec++;
				}
			}
		}
		return mec;
	}

	public double getIntraClusterDistance(List<RefCpG> twoClusterRefCpGList) {
		double distance = 0;
		for (MappedRead mappedRead : mappedReadList) {
			double readAvgDistance = 0, validCpGCount = 0;
			for (CpG cpG : mappedRead.getCpgList()) {
				for (RefCpG refCpG : twoClusterRefCpGList) {
					if (cpG.getPos() == refCpG.getPos()) {
						double majorMethylLevel = refCpGMap.get(cpG.getPos()).getMethylLevel();
						if (cpG.getMethylStatus() == MethylStatus.C) {
							readAvgDistance += 1 - majorMethylLevel;
							validCpGCount++;
						} else if (cpG.getMethylStatus() == MethylStatus.T) {
							readAvgDistance += majorMethylLevel;
							validCpGCount++;
						}
					}
				}
			}
			readAvgDistance /= validCpGCount;
			distance += readAvgDistance;
		}
		return distance / mappedReadList.size();
	}

	public int getCpGSum() {
		int cpgSum = 0;
		for (MappedRead mappedRead : mappedReadList) {
			cpgSum += mappedRead.getCpgList().size();
		}
		return cpgSum;
	}
}
