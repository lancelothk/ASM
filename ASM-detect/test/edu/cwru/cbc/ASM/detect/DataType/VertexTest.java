package edu.cwru.cbc.ASM.detect.DataType;

import edu.cwru.cbc.ASM.commons.Methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static org.testng.AssertJUnit.assertEquals;


public class VertexTest {
	public MappedRead mappedRead;
	public Vertex vertex;

	@BeforeMethod
	public void setUp() throws Exception {
		mappedRead = new MappedRead("test", '+', 0, "ACGTGTGCAG", "test-1");
		mappedRead.generateCpGsInRead(MethylationUtils.extractCpGSite("ACGCGTGCAG", 0).stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG)));
		vertex = new Vertex(mappedRead);
	}

	@Test
	public void testGetMECScore() throws Exception {
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), vertex.getRefCpGMap().get(1), 4);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), vertex.getRefCpGMap().get(1), 8);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), vertex.getRefCpGMap().get(3), 5);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), vertex.getRefCpGMap().get(3), 8);
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