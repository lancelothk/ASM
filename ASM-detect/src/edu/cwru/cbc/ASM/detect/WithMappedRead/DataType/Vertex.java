package edu.cwru.cbc.ASM.detect.WithMappedRead.DataType;

import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MethylStatus;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.util.*;

/**
 * Created by lancelothk on 5/14/14.
 * Vertex in graph.
 */
public class Vertex {
	private String id;
	private List<Edge> adjEdges;
	private List<MappedRead> mappedReadList;
	private Map<Integer, RefCpG> refCpGMap; // position-RefCpG
	private int methylPolarity;

	public Vertex(MappedRead mappedRead) {
		this.id = mappedRead.getId();
		this.adjEdges = new ArrayList<>();
		this.mappedReadList = new ArrayList<>();
		this.mappedReadList.add(mappedRead);
		this.refCpGMap = new HashMap<>();
		this.methylPolarity = mappedRead.getMethylPolarity();
		this.addCpG(mappedRead.getCpgList());
	}

	public void addCpG(Collection<CpG> cpgList) {
		cpgList.stream().filter(cpg -> !refCpGMap.containsKey(cpg.getPos())).forEach(
				cpg -> refCpGMap.put(cpg.getPos(), new RefCpG(cpg.getPos(), cpg.getMethylStatus())));
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

	public int getFirstCpGPos() {
		return refCpGMap.values().stream().min((c1, c2) -> c1.getPos() - c2.getPos()).get().getPos();
	}

	public void addEdge(Edge edge) {
		this.adjEdges.add(edge);
	}

	public void removeEdge(Edge edge) {
		this.adjEdges.remove(edge);
	}

	public void addMappedRead(List<MappedRead> mappedReadList) {
		mappedReadList.forEach((read) -> methylPolarity += read.getMethylPolarity());
		this.mappedReadList.addAll(mappedReadList);
	}

	public String getId() {
		return id;
	}

	public List<Edge> getAdjEdges() {
		return adjEdges;
	}

	public List<MappedRead> getMappedReadList() {
		return mappedReadList;
	}

	public int getMethylPolarity() {
		return methylPolarity;
	}

	public double getMECScore() {
		double mec = 0;
		for (MappedRead mappedRead : mappedReadList) {
			for (CpG cpG : mappedRead.getCpgList()) {
					MethylStatus refMethylStatus = refCpGMap.get(cpG.getPos()).getMajorMethylStatus();
				// TODO double check the N in calculation of MEC
					if (refMethylStatus == MethylStatus.N) {
						mec += 0.5;
					} else if (cpG.getMethylStatus() != refMethylStatus && cpG.getMethylStatus() != MethylStatus.N) {
						mec++;
					}
			}
		}
		return mec;
	}

	public int getCpGSum() {
		int cpgSum = 0;
		for (MappedRead mappedRead : mappedReadList) {
			cpgSum += mappedRead.getCpgList().size();
		}
		return cpgSum;
	}
}
