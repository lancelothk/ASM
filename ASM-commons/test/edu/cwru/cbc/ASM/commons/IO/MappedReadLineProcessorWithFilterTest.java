package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import org.testng.annotations.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertTrue;

public class MappedReadLineProcessorWithFilterTest {

	@Test
	public void testProcessLine() throws Exception {
		String refString = "CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG";
		List<RefCpG> refCpGList = extractCpGSite(refString, 17207805);

		// mapped read is 0based start.
		String mappedReadStr1 = "20\t+\t17207806\t17207874\tGTATTTCGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		String mappedReadStr3 = "20\t+\t17207806\t17207874\tGTATTTNGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815506";
		// contains only one CpG, shouldn't be added to mappredReadList
		String mappedReadStr2 = "20\t-\t17207806\t17207814\tACATTCGCA\t3458074"; // not the original sequence
		MappedReadLineProcessorWithFilter mp = new MappedReadLineProcessorWithFilter(refCpGList, 2, refString.length());

		mp.processLine(mappedReadStr1);
		mp.processLine(mappedReadStr2);
		mp.processLine(mappedReadStr3);

		// check if mappedReadList exclude read with only one cpg
		assertEquals("incorrect mappedReadList size!", 2, mp.mappedReadList.size());
		assertTrue("mappredReadList shouldn't contains 1 CpG read!", mp.mappedReadList.stream().noneMatch(r -> r.getId().equals("3458074")));

		// check if N CpG included in refMap and mappedRead
		assertEquals("incorrect refMap CpG number!", 1, mp.refMap.get(17207812).getCpGList().size());
		assertEquals("incorrect refMap CpG number!", 2, mp.refMap.get(17207838).getCpGList().size());
		assertEquals("incorrect mappedRead CpG number!", 3, mp.mappedReadList.get(0).getCpgList().size());
		assertEquals("incorrect mappedRead CpG number!", 2, mp.mappedReadList.get(1).getCpgList().size());
	}
}