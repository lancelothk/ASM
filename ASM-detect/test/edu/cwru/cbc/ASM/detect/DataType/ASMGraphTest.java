package edu.cwru.cbc.ASM.detect.DataType;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;
import edu.cwru.cbc.ASM.commons.Read.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.ReflectionUtils;
import edu.cwru.cbc.ASM.detect.DetectionUtils;
import org.junit.Test;

import java.io.File;
import java.util.List;
import java.util.Map;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;
import static org.junit.Assert.assertEquals;

public class ASMGraphTest {

	@Test
	public void testCluster() throws Exception {
		String homeDirectory = System.getProperty("user.home");
		File inputFile = new File(homeDirectory + "/IdeaProjects/ASM/ASM-detect/testData/chr20-56859339-56859404");

		int startPos = Integer.parseInt(inputFile.getName().split("-")[1]);

		// load input
		String reference = DetectionUtils.readRefFromIntervalFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor(refCpGList));
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