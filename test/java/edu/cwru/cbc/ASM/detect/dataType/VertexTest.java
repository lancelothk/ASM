package edu.cwru.cbc.ASM.detect.dataType;

import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.lang3.reflect.FieldUtils;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static org.testng.AssertJUnit.assertEquals;


public class VertexTest {
	public MappedRead mappedRead;
	public Vertex vertex;

	@BeforeMethod
	public void setUp() throws Exception {
		mappedRead = new MappedRead("test", '+', 0, "ACGTGTGCAG", "1");
		List<RefCpG> refCpGList = MethylationUtils.extractCpGSite("ACGCGTGCAG", 0);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		mappedRead.generateCpGsInRead(refMap);
		vertex = new Vertex(mappedRead);
	}

	@Test
	public void testGetMECScore() throws Exception {
		FieldUtils.writeField(vertex.getRefCpGMap().get(1), "methylCount", 4, true);
		FieldUtils.writeField(vertex.getRefCpGMap().get(1), "coveredCount", 8, true);
		FieldUtils.writeField(vertex.getRefCpGMap().get(3), "methylCount", 5, true);
		FieldUtils.writeField(vertex.getRefCpGMap().get(3), "coveredCount", 8, true);
		assertEquals("incorrect MEC score", 1.5, vertex.getMECScore(), 0.00001);
	}

	@Test
	public void testAddCpG() throws Exception {
		assertEquals("incorrect RefCpG size", vertex.getRefCpGMap().size(), 2);
	}

	@Test
	public void testAddRefCpG() throws Exception {
		List<RefCpG> refCpGList = new ArrayList<>();
		refCpGList.add(new RefCpG(1, MethylStatus.T));
		refCpGList.add(new RefCpG(3, MethylStatus.T));
		refCpGList.add(new RefCpG(5, MethylStatus.C));
		vertex.addRefCpG(refCpGList);
		assertEquals("incorrect RefCpG size", vertex.getRefCpGMap().size(), 3);
	}
}