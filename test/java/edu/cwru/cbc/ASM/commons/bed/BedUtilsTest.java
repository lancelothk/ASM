package edu.cwru.cbc.ASM.commons.bed;

import edu.cwru.cbc.ASM.commons.GenomicInterval.BedInterval;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Map;

import static org.testng.AssertJUnit.assertEquals;


/**
 * Created by lancelothk on 4/2/15.
 * Test Bed utils
 */
public class BedUtilsTest {

	@Test
	public void testReadBedRegions() throws Exception {
		Map<String, List<BedInterval>> bedRegionsWithLabel = BedUtils.readBedRegions(
				getClass().getClassLoader().getResource("bedWithLabel.bed").getFile(), true);

		assertEquals("Bed region size incorrect!", 11, bedRegionsWithLabel.get("chr20").size());
		assertEquals("incorrect label!", true, bedRegionsWithLabel.get("chr20").get(0).isPositive());


		// TODO add error cases.
	}

}