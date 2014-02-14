package ASM.Execution;

import ASM.DataType.CpGIsland;
import org.junit.Test;

import java.util.List;

public class MatchCpGIslandTest {

	@Test
	public void testReadCpGIslands() throws Exception {
		List<CpGIsland> cpgIslandList = MatchCpGIsland.readCpGIslands("/media/ke/win-data/Dataset/cpgIslandExt_hg18_UCSCGB_chr6.txt");
		for (int i = 0; i < 5; i++) {
			System.out.println(cpgIslandList.get(i).toString());
		}
	}

	@Test
	public void testMatchMappedReadToCpGIsland() throws Exception {
		MatchCpGIsland.matchMappedReadToCpGIsland("/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/reads_bs_i90_r1.mapped_chr6_head100000",
				"/media/ke/win-data/Dataset/cpgIslandExt_hg18_UCSCGB_chr6.txt", "/home/ke/test");
	}

}
