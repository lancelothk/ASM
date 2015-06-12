package edu.cwru.cbc.ASM.commons.IO;

import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.Sequence.MappedRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils.extractCpGSite;
import static org.testng.AssertJUnit.assertEquals;

public class MappedReadLineProcessorTest {

	@DataProvider(name = "invalidCharacters")
	public static Object[][] invalidCharacters() {
		String invalidCharacters = "acgtnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm@#$%+_=-~!^&*()/\\?><,:;\"\' ";
		Object[][] result = new Object[invalidCharacters.length()][];
		for (int i = 0; i < invalidCharacters.length(); i++) {
			result[i] = new Object[]{invalidCharacters.charAt(i)};
		}
		return result;
	}

	@DataProvider(name = "validCharacters")
	public static Object[][] validCharacters() {
		String validCharacters = "ACGTN.";
		Object[][] result = new Object[validCharacters.length()][];
		for (int i = 0; i < validCharacters.length(); i++) {
			result[i] = new Object[]{validCharacters.charAt(i)};
		}
		return result;
	}

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

	@Test(dataProvider = "validCharacters")
	public void test_processLine_validCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(new ArrayList<>());
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
		List<MappedRead> mappedReadList = mlp.getResult();
		assertEquals("incorrect character!", c, mappedReadList.get(0).getSequence().charAt(0));
	}

	@Test(expectedExceptions = RuntimeException.class, dataProvider = "invalidCharacters")
	public void test_processLine_invalidCharacter(final char c) throws Exception {
		MappedReadLineProcessor mlp = new MappedReadLineProcessor(new ArrayList<>());
		String mappedReadStr1 = "20\t+\t17207806\t17207874\t" + c + "\t815505";
		mlp.processLine(mappedReadStr1);
	}
}