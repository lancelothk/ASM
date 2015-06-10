package edu.cwru.cbc.ASM.detect.DataType;

import edu.cwru.cbc.ASM.commons.Methylation.CpG;
import edu.cwru.cbc.ASM.commons.Methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.junit.Before;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class VertexTest {
	public MappedRead mappedRead;
	public Vertex vertex;

	@Before
	public void setUp() throws Exception {
		mappedRead = new MappedRead("test", '+', 0, "ACGTGTGCAG", "test-1");
		mappedRead.addCpG(new CpG(mappedRead, new RefCpG(1), MethylStatus.C));
		mappedRead.addCpG(new CpG(mappedRead, new RefCpG(3), MethylStatus.T));
		vertex = new Vertex(mappedRead);
	}

	@org.junit.Test
	public void testGetMECScore() throws Exception {
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), vertex.getRefCpGMap().get(1), 4);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), vertex.getRefCpGMap().get(1), 8);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("methylCount"), vertex.getRefCpGMap().get(3), 5);
		ReflectionUtils.setPrivateField(RefCpG.class.getDeclaredField("coveredCount"), vertex.getRefCpGMap().get(3), 8);
		assertEquals("incorrect MEC score", 1.5, vertex.getMECScore(), 0.00001);
	}

	@org.junit.Test
	public void testAddCpG() throws Exception {
		assertEquals("incorrect RefCpG size", vertex.getRefCpGMap().size(), 2);
	}

	@org.junit.Test
	public void testAddRefCpG() throws Exception {
		List<RefCpG> refCpGList = new ArrayList<>();
		refCpGList.add(new RefCpG(1, MethylStatus.T));
		refCpGList.add(new RefCpG(3, MethylStatus.T));
		refCpGList.add(new RefCpG(5, MethylStatus.C));
		vertex.addRefCpG(refCpGList);
		assertEquals("incorrect RefCpG size", vertex.getRefCpGMap().size(), 3);
	}
}