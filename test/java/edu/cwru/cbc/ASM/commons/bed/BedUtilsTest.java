package edu.cwru.cbc.ASM.commons.bed;

import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import org.testng.annotations.Test;

import java.net.URL;
import java.util.List;
import java.util.Map;

import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertNotNull;


/**
 * Created by lancelothk on 4/2/15.
 * Test Bed utils
 */
public class BedUtilsTest {

	@Test
	public void testReadBedRegions() throws Exception {
		URL file = getClass().getClassLoader().getResource("bedWithLabel.bed");
		assertNotNull(file);
		Map<String, List<BedInterval>> bedRegionsWithLabel = BedUtils.readBedRegions(file.getFile(), true);

		assertEquals("Bed region size incorrect!", 11, bedRegionsWithLabel.get("chr20").size());
		assertEquals("incorrect label!", true, bedRegionsWithLabel.get("chr20").get(0).isPositive());


		// TODO add error cases.
	}

}