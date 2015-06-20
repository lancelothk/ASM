package edu.cwru.cbc.ASM.commons.sequence;

import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.stream.Collectors;

import static org.testng.Assert.assertEquals;

public class MappedReadTest {
	private Map<Integer, RefCpG> refCpGMap;

	@BeforeMethod
	public void setUp() throws Exception {
		refCpGMap = MethylationUtils.extractCpGSite("CGATCGACGACG", 0).stream().collect(
				Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
	}

	@Test
	public void testGenerateCpGsInRead() throws Exception {
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", 1);
		MappedRead minusStrandRead = new MappedRead("minus", '-', 0, "CTATGCTCNTGC", 2);
		MappedRead plusStrandRead_partialCpG = new MappedRead("plus", '+', 1, "GATCGANGAT", 3);
		MappedRead minusStrandRead_partialCpG = new MappedRead("minus", '-', 1, "CTAGTTGNTG", 4);

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
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", 1);
		MappedRead minusStrandRead = new MappedRead("minus", '-', 0, "CTATGCTCNTGC", 2);
		assertEquals(plusStrandRead.getComplementarySequence(), "ACTAGCTNCTGC");
		assertEquals(minusStrandRead.getComplementarySequence(), "GATACGAGNACG");
	}

	@Test
	public void testToVisualizationString() throws Exception {
		MappedRead read = new MappedRead("6", '+', 10, "1234567890", 1);
		assertEquals(read.toVisualizationString(5), "6\t+\t10\t19\t.....1234567890\t1");
	}

	@Test
	public void testToSimulationString() throws Exception {
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", 1);
		plusStrandRead.generateCpGsInRead(refCpGMap);
		plusStrandRead.getCpgList().get(0).setMethylStatus(MethylStatus.C);
		plusStrandRead.getCpgList().get(1).setMethylStatus(MethylStatus.T);
		plusStrandRead.getCpgList().get(2).setMethylStatus(MethylStatus.N);
		assertEquals(plusStrandRead.toSimulationString(), "plus\t+\t0\t11\tCGATTGANGACG\t1\n");

		MappedRead minusStrandRead = new MappedRead("minus", '-', 0, "CTATGCTCNTGC", 2);
		minusStrandRead.generateCpGsInRead(refCpGMap);
		minusStrandRead.getCpgList().get(0).setMethylStatus(MethylStatus.C);
		minusStrandRead.getCpgList().get(1).setMethylStatus(MethylStatus.C);
		minusStrandRead.getCpgList().get(2).setMethylStatus(MethylStatus.T);
		assertEquals(minusStrandRead.toSimulationString(), "minus\t-\t0\t11\tGCATGCTCNTGT\t2\n");
	}

	@Test
	public void testToMRFormatString() throws Exception {
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGAT", 1);
		assertEquals(plusStrandRead.toMRFormatString(3, 'a'), "plus\t0\t3\t1\t3\t+\tTGAT\taaaa\n");
	}

	@Test
	public void testMappedReadString() throws Exception {
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", 1);
		assertEquals(plusStrandRead.toMappedReadString(), "plus\t+\t0\t11\tTGATCGANGACG\t1\n");
	}

	@Test
	public void testMappedReadString_withBound() throws Exception {
		MappedRead plusStrandRead = new MappedRead("plus", '+', 0, "TGATCGANGACG", 1);
		assertEquals(plusStrandRead.toMappedReadString(1, 4), "plus\t+\t1\t4\tGATC\t1\n");
	}
}