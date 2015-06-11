package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import org.testng.annotations.Test;

import java.util.List;

import static edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;

public class MappedReadLineProcessorTest {

	@Test
	public void testProcessLine() throws Exception {
		String refString = "CGCATTCCGATGCAGAATGTCCTTCATGAGAGGCGACTTTTTAGGACTTTTAATCTGCGTTCAAATCAATTAATAGTTTGACG";
		List<RefCpG> refCpGList = extractCpGSite(refString, 17207805);

		// mapped read is 0based start.
		String mappedReadStr1 = "20\t+\t17207806\t17207874\tGTATTTCGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815505";
		// read with one CpG
		String mappedReadStr2 = "20\t-\t17207806\t17207814\tACATTCGCA\t3458074";// not the original sequence
		String mappedReadStr3 = "20\t+\t17207806\t17207874\tGTATTTNGATGTAGAATGTTTTTTATGAGAGGTGATTTTTTAGGATTTTTAATTTGTGTTTAAATTAAT\t815506";

		MappedReadLineProcessor mp = new MappedReadLineProcessor(refCpGList);
		mp.processLine(mappedReadStr1);
		mp.processLine(mappedReadStr2);
		mp.processLine(mappedReadStr3);

		// check if mappedReadList include all reads
		assertEquals("incorrect mappedReadList size!", 3, mp.mappedReadList.size());

		// check if N CpG included in refMap and mappedRead
		assertEquals("incorrect refMap CpG number!", 2, mp.refMap.get(17207812).getCpGList().size());
		assertEquals("incorrect refMap CpG number!", 2, mp.refMap.get(17207838).getCpGList().size());
		// here the order of read is indirectly tested.
		assertEquals("incorrect mappedRead CpG number!", 3, mp.mappedReadList.get(0).getCpgList().size());
		assertEquals("incorrect mappedRead CpG number!", 1, mp.mappedReadList.get(1).getCpgList().size());
		assertEquals("incorrect mappedRead CpG number!", 2, mp.mappedReadList.get(2).getCpgList().size());
	}
}