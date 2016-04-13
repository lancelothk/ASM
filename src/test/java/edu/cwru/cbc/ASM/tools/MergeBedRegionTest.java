package edu.cwru.cbc.ASM.tools;

import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.testng.Assert.assertEquals;

/**
 * Created by lancelothk on 7/6/15.
 */
public class MergeBedRegionTest {
	Map<String, List<BedInterval>> bedRegionMap;

	@BeforeMethod
	public void setUp() throws Exception {
		bedRegionMap = new HashMap<>();
		bedRegionMap.put("chr1", new ArrayList<>());
		bedRegionMap.put("chr2", new ArrayList<>());
		bedRegionMap.put("chr3", new ArrayList<>());
		// continuous regions with in range 30
		bedRegionMap.get("chr1").add(new BedInterval("chr1", 20, 50, "test11"));
		bedRegionMap.get("chr1").add(new BedInterval("chr1", 80, 100, "test12"));
		bedRegionMap.get("chr1").add(new BedInterval("chr1", 120, 150, "test13"));
		// region gap > 30 and regions with overlap
		bedRegionMap.get("chr2").add(new BedInterval("chr2", 20, 50, "test21"));
		bedRegionMap.get("chr2").add(new BedInterval("chr2", 100, 150, "test22"));
		bedRegionMap.get("chr2").add(new BedInterval("chr2", 120, 180, "test23"));
		// region inside another region
		bedRegionMap.get("chr3").add(new BedInterval("chr3", 100, 150, "test31"));
		bedRegionMap.get("chr3").add(new BedInterval("chr3", 120, 140, "test32"));
	}

	@Test
	public void test_MergeBedRegion() throws Exception {
		Method mergeBedRegions = MergeBedRegion.class.getDeclaredMethod("mergeBedRegions", Map.class, int.class);
		mergeBedRegions.setAccessible(true);
		mergeBedRegions.invoke(null, bedRegionMap, 30);
		assertEquals(bedRegionMap.get("chr1").size(), 1);
		assertEquals(bedRegionMap.get("chr2").size(), 2);
		assertEquals(bedRegionMap.get("chr3").size(), 1);
	}
}