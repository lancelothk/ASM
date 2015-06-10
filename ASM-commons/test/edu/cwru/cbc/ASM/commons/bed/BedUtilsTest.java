package edu.cwru.cbc.ASM.commons.bed;

import edu.cwru.cbc.ASM.commons.GenomicInterval.BedInterval;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Created by lancelothk on 4/2/15.
 * Test Bed utils
 */
public class BedUtilsTest {

	@Test
	public void testReadBedRegions() throws Exception {
		String bedWithLabel = "testData/bedWithLabel.bed";

		Map<String, List<BedInterval>> bedRegionsWithLabel = BedUtils.readBedRegions(bedWithLabel, true);

		assertEquals("Bed region size incorrect!", 11, bedRegionsWithLabel.get("chr20").size());
		assertEquals("incorrect label!", true, bedRegionsWithLabel.get("chr20").get(0).isPositive());


		// TODO add error cases.
	}

}