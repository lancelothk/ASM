package edu.cwru.cbc.ASM.detect.DataType;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.IO.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import edu.cwru.cbc.ASM.detect.DetectionUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;

public class ASMGraphTest {

	@Test
	public void testCluster() throws Exception {
		File inputFile = new File(getClass().getClassLoader().getResource("chr20-56859339-56859404").getFile());

		int startPos = Integer.parseInt(inputFile.getName().split("-")[1]);

		// load input
		String reference = DetectionUtils.readRefFromIntervalFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor());
		Map<Integer, RefCpG> refMap = refCpGList.stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refMap));

		ASMGraph asmGraph = new ASMGraph(mappedReadList);
		List edgeList = (List) (ReflectionUtils.getPrivateField(ASMGraph.class.getDeclaredField("edgeList"), asmGraph));
		Map vertexMap = (Map) (ReflectionUtils.getPrivateField(ASMGraph.class.getDeclaredField("vertexMap"), asmGraph));
		assertEquals("incorrect edge number", 43, edgeList.size());
		assertEquals("incorrect edge number", 10, vertexMap.size());

		for (Object o : edgeList) {
			Edge e = (Edge) o;
			if ((e.getLeft().getId().equals("10040340") && e.getRight().getId().equals("5087521")) ||
					(e.getLeft().getId().equals("5087521") && e.getRight().getId().equals("10040340"))) {
				assertEquals("incorrect edge weight!", -2, e.getWeight(), 0.0001);
			}
		}

		asmGraph.cluster();

		assertEquals("incorrect cluster number!", 2, asmGraph.getClusterResult().size());
		assertEquals("incorrect cluster number!", 1, asmGraph.getMECSum(), 0.0001);
	}
}