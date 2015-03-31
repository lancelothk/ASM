package edu.cwru.cbc.ASM.commons.DataType;

import org.junit.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

public class MappedReadLineProcessorTest {

	@Test
	public void testProcessRead() throws Exception {
		List<RefCpG> refCpGList = extractCpGSite(
				"CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG", 17207805);


		String mappedReadStr1 = "20\t+\t17207806\t17207874\tGTATTTCGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		String mappedReadStr3 = "20\t+\t17207806\t17207874\tGTATTTNGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		// contains only one CpG, shouldn't be added to mappredReadList
		String mappedReadStr2 = "20\t-\t17207806\t17207814\tACATTCCGA\t3458074";
		MappedReadLineProcessor mp = new MappedReadLineProcessor(refCpGList, 2);

		MappedRead mappedRead1 = mp.processRead(mappedReadStr1);
		MappedRead mappedRead2 = mp.processRead(mappedReadStr2);
		MappedRead mappedRead3 = mp.processRead(mappedReadStr3);

		// check if mappedReadList exclude read with only one cpg
		assertEquals("incorrect mappedReadList size!", mp.mappedReadList.size(), 2);
		assertFalse("mappredReadList shouldn't contains 1 CpG read!", mp.mappedReadList.contains(mappedRead2));

		// check if N CpG included in refMap and mappedRead
		assertEquals("incorrect refMap CpG number!", mp.refMap.get(17207812).getCpGList().size(), 1);
		assertEquals("incorrect refMap CpG number!", mp.refMap.get(17207838).getCpGList().size(), 2);
		assertEquals("incorrect mappedRead CpG number!", mappedRead1.getCpgList().size(), 3);
		assertEquals("incorrect mappedRead CpG number!", mappedRead3.getCpgList().size(), 2);
	}
}