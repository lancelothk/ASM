package edu.cwru.cbc.ASM.commons.Sequence;

import edu.cwru.cbc.ASM.commons.Methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.Methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.Methylation.RefCpG;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.stream.Collectors;

import static org.testng.Assert.assertEquals;

public class MappedReadTest {
	private MappedRead plusStrandRead;
	private MappedRead minusStrandRead;
	private MappedRead plusStrandRead_partialCpG;
	private MappedRead minusStrandRead_partialCpG;
	private Map<Integer, RefCpG> refCpGMap;

	@BeforeMethod
	public void setUp() throws Exception {
		refCpGMap = MethylationUtils.extractCpGSite("CGATCGACGACG", 0).stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
		plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", "read_plus1");
		minusStrandRead = new MappedRead("minus", '-', 0, "CTATGCTCNTGC", "read_minus1");
		plusStrandRead_partialCpG = new MappedRead("plus", '+', 1, "GATCGANGAT", "read_plus2");
		minusStrandRead_partialCpG = new MappedRead("minus", '-', 1, "CTAGTTGNTG", "read_minus2");
	}

	@Test
	public void testGenerateCpGsInRead() throws Exception {
		plusStrandRead.generateCpGsInRead(refCpGMap);
		minusStrandRead.generateCpGsInRead(refCpGMap);
		plusStrandRead_partialCpG.generateCpGsInRead(refCpGMap);
		minusStrandRead_partialCpG.generateCpGsInRead(refCpGMap);

		assertEquals(3, plusStrandRead.getCpgList().size());
		assertEquals(MethylStatus.T, plusStrandRead.getCpgList().get(0).getMethylStatus());
		assertEquals(MethylStatus.C, plusStrandRead.getCpgList().get(1).getMethylStatus());
		assertEquals(MethylStatus.C, plusStrandRead.getCpgList().get(2).getMethylStatus());

		assertEquals(3, minusStrandRead.getCpgList().size());
		assertEquals(MethylStatus.T, minusStrandRead.getCpgList().get(0).getMethylStatus());
		assertEquals(MethylStatus.C, minusStrandRead.getCpgList().get(1).getMethylStatus());
		assertEquals(MethylStatus.C, minusStrandRead.getCpgList().get(2).getMethylStatus());

		assertEquals(2, plusStrandRead_partialCpG.getCpgList().size());
		assertEquals(MethylStatus.C, plusStrandRead_partialCpG.getCpgList().get(0).getMethylStatus());
		assertEquals(MethylStatus.T, plusStrandRead_partialCpG.getCpgList().get(1).getMethylStatus());

		assertEquals(2, minusStrandRead_partialCpG.getCpgList().size());
		assertEquals(MethylStatus.C, minusStrandRead_partialCpG.getCpgList().get(0).getMethylStatus());
		assertEquals(MethylStatus.T, minusStrandRead_partialCpG.getCpgList().get(1).getMethylStatus());
	}

	@Test
	public void testGetComplementarySequence() throws Exception {
		assertEquals(plusStrandRead.getComplementarySequence(), "ACTAGCTNCTGC");
		assertEquals(minusStrandRead.getComplementarySequence(), "GATACGAGNACG");
	}

	@Test
	public void testToString() throws Exception {
		MappedRead read = new MappedRead("6", '+', 10, "1234567890", "test");
		assertEquals(read.toVisualizationString(5), "6\t+\t10\t20\t.....1234567890\ttest");
	}
}