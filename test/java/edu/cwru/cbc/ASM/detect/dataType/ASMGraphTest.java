package edu.cwru.cbc.ASM.detect.dataType;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.io.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.lang3.reflect.FieldUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.net.URL;
import java.util.List;
import java.util.Map;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertNotNull;

public class ASMGraphTest {

	@Test
	public void testCluster() throws Exception {
		URL file = getClass().getClassLoader().getResource("chr20-56859339-56859404");
		assertNotNull(file);
		File inputFile = new File(file.getFile());

		int startPos = Integer.parseInt(inputFile.getName().split("-")[1]);

		// load input
		String reference = IOUtils.readRefFromIntervalReadsFile(inputFile);
		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor());
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refMap));

		ASMGraph asmGraph = new ASMGraph(mappedReadList);
		List edgeList = (List) FieldUtils.readDeclaredField(asmGraph, "edgeList", true);
		Map vertexMap = (Map) FieldUtils.readDeclaredField(asmGraph, "vertexMap", true);
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