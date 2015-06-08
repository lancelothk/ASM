package edu.cwru.cbc.ASM.commons.Read;

import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import org.junit.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.CommonsUtils.extractCpGSite;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

public class MappedReadLineProcessorTest {

	@Test
	public void testProcessRead() throws Exception {
		String refString = "CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG";
		List<RefCpG> refCpGList = extractCpGSite(refString, 17207805);


		String mappedReadStr1 = "20\t+\t17207806\t17207874\tGTATTTCGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		String mappedReadStr3 = "20\t+\t17207806\t17207874\tGTATTTNGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		// contains only one CpG, shouldn't be added to mappredReadList
		String mappedReadStr2 = "20\t-\t17207806\t17207814\tACATTCCGA\t3458074";
		MappedReadLineProcessorWithSummary mp = new MappedReadLineProcessorWithSummary(refCpGList, 2, refString.length());

		MappedRead mappedRead1 = mp.processRead(mappedReadStr1);
		MappedRead mappedRead2 = mp.processRead(mappedReadStr2);
		MappedRead mappedRead3 = mp.processRead(mappedReadStr3);
		mp.updateRefCpG(mappedRead1);
		mp.updateRefCpG(mappedRead2);
		mp.updateRefCpG(mappedRead3);

		// check if mappedReadList exclude read with only one cpg
		assertEquals("incorrect mappedReadList size!", 2, mp.mappedReadList.size());
		assertFalse("mappredReadList shouldn't contains 1 CpG read!", mp.mappedReadList.contains(mappedRead2));

		// check if N CpG included in refMap and mappedRead
		assertEquals("incorrect refMap CpG number!", 1, mp.refMap.get(17207812).getCpGList().size());
		assertEquals("incorrect refMap CpG number!", 2, mp.refMap.get(17207838).getCpGList().size());
		assertEquals("incorrect mappedRead CpG number!", 3, mappedRead1.getCpgList().size());
		assertEquals("incorrect mappedRead CpG number!", 2, mappedRead3.getCpgList().size());
	}
}