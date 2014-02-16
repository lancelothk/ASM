package ASM.Utils;

import ASM.DataType.CpGIsland;
import org.junit.Test;

import java.util.List;

/**
 * Created by lancelothk on 2/16/14.
 */
public class ReadTextFileTest {

	@Test
	public void testReadTextFileWithDelimiter() throws Exception {
		List<CpGIsland> cpGIslandList = ReadTextFile.ReadTextFileWithDelimiter(CpGIsland.class.getDeclaredConstructor(String.class, long.class, long.class, String.class, int.class, int.class, int.class, double.class, double.class,double.class)
		,"/home/lancelothk/Documents/cpgIslandExt_hg18_UCSCGB_chr6.txt", "\t", false);
		for (CpGIsland cpGIsland : cpGIslandList) {
			System.out.println(cpGIsland.toString());
		}
	}
}