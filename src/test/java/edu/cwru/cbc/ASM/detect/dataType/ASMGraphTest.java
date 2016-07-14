package edu.cwru.cbc.ASM.detect.dataType;

public class ASMGraphTest {

//	@Test
//	public void testCluster() throws Exception {
//		URL file = getClass().getClassLoader().getResource("chr20-56859339-56859404");
//		AssertJUnit.assertNotNull(file);
//		File inputFile = new File(file.getFile());
//		int startPos = Integer.parseInt(inputFile.getName().split("-")[1]);
//
//		// load input
//		String reference = IOUtils.readRefFromIntervalReadsFile(inputFile);
//		List<RefCpG> refCpGList = extractCpGSite(reference, startPos);
//		List<MappedRead> mappedReadList = Files.asCharSource(inputFile, Charsets.UTF_8)
//				.readLines(new MappedReadLineProcessor());
//		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
//		for (RefCpG refCpG : refCpGList) {
//			refMap.put(refCpG.getPos(), refCpG);
//		}
//		mappedReadList.forEach(mr -> mr.generateCpGsInRead(refMap));
//
//		ASMGraph asmGraph = new ASMGraph(mappedReadList);
//		Set edgeList = (Set) FieldUtils.readDeclaredField(asmGraph, "edgeSet", true);
//		Map vertexMap = (Map) FieldUtils.readDeclaredField(asmGraph, "vertexMap", true);
//		AssertJUnit.assertEquals("incorrect edge number", 43, edgeList.size());
//		AssertJUnit.assertEquals("incorrect edge number", 10, vertexMap.size());
//
//		for (Object o : edgeList) {
//			Edge e = (Edge) o;
//			if ((e.getLeft().getId().equals("10040340") && e.getRight().getId().equals("5087521")) ||
//					(e.getLeft().getId().equals("5087521") && e.getRight().getId().equals("10040340"))) {
//				AssertJUnit.assertEquals("incorrect edge weight!", -2, e.getWeight(), 0.0001);
//			}
//		}
//
//		asmGraph.cluster();
//
//		AssertJUnit.assertEquals("incorrect cluster number!", 2, asmGraph.getClusterResult().size());
//		AssertJUnit.assertEquals("incorrect cluster number!", 1, asmGraph.getMECSum(), 0.0001);
//	}
}