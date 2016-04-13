package edu.cwru.cbc.ASM.commons.bed;

import edu.cwru.cbc.ASM.commons.genomicInterval.BedInterval;
import org.testng.annotations.Test;

import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import static org.testng.AssertJUnit.*;


/**
 * Created by lancelothk on 4/2/15.
 * Test Bed utils
 */
public class BedUtilsTest {

	@Test
	public void testReadBedRegions() throws Exception {
		URL file = getClass().getClassLoader().getResource("multiChrBedFile.bed");
		assertNotNull(file);
		Map<String, List<BedInterval>> bedRegionsMap = BedUtils.readBedRegions(file.getFile());
		assertEquals("Bed region size incorrect!", 2, bedRegionsMap.get("chr1").size());
		assertEquals("Bed region size incorrect!", 2, bedRegionsMap.get("chr3").size());
		assertEquals("Bed region size incorrect!", 1, bedRegionsMap.get("chr6").size());
	}

	@Test
	public void testReadBedRegions_3col() throws Exception {
		URL file = getClass().getClassLoader().getResource("multiChrBedFile_3col.bed");
		assertNotNull(file);
		Map<String, List<BedInterval>> bedRegionsMap = BedUtils.readBedRegions(file.getFile());
		assertEquals("Bed region size incorrect!", 2, bedRegionsMap.get("chr1").size());
		assertEquals("Bed region size incorrect!", 2, bedRegionsMap.get("chr3").size());
		assertEquals("Bed region size incorrect!", 1, bedRegionsMap.get("chr6").size());
	}

	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*less than 3 columns in bed file.*")
	public void testReadBedRegions_2col() throws Exception {
		URL file = getClass().getClassLoader().getResource("invalid_2col.bed");
		assertNotNull(file);
		BedUtils.readBedRegions(file.getFile());
	}


	@Test(expectedExceptions = RuntimeException.class,
			expectedExceptionsMessageRegExp = ".*Bed file contains regions from multiple chromosomes.*")
	public void testReadSingleChrBedRegions() throws Exception {
		URL file = getClass().getClassLoader().getResource("singleChrBedFile.bed");
		assertNotNull(file);
		List<BedInterval> bedRegions = BedUtils.readSingleChromBedRegions(file.getFile());
		assertEquals("Bed region size incorrect!", 11, bedRegions.size());

		file = getClass().getClassLoader().getResource("multiChrBedFile.bed");
		assertNotNull(file);
		BedUtils.readSingleChromBedRegions(file.getFile());
	}


	@Test
	public void testIntersect() throws Exception {
		BedInterval a1 = new BedInterval("chr1", 5, 20, "1");
		BedInterval b1 = new BedInterval("chr1", 0, 9, "1");
		BedInterval b2 = new BedInterval("chr1", 10, 14, "2");
		BedInterval b3 = new BedInterval("chr1", 15, 24, "3");
		List<BedInterval> a = new ArrayList<>();
		a.add(a1);
		List<BedInterval> b = new ArrayList<>();
		b.add(b1);
		b.add(b2);
		b.add(b3);

		BedInterval a1_2 = new BedInterval("chr2", 5, 20, "1");
		BedInterval b1_2 = new BedInterval("chr2", 0, 9, "1");
		BedInterval b2_2 = new BedInterval("chr2", 10, 14, "2");
		BedInterval b3_2 = new BedInterval("chr2", 15, 24, "3");
		a.add(a1_2);
		b.add(b1_2);
		b.add(b2_2);
		b.add(b3_2);

		Collection<BedInterval> intersections = BedUtils.intersect(a, b);
		assertEquals(6, intersections.size());

		ArrayList<BedInterval> intersectionList = new ArrayList<>(intersections);
		intersectionList.sort(BedInterval::compareTo);
		assertEquals("chr1", intersectionList.get(0).getChr());
		assertEquals("chr1", intersectionList.get(1).getChr());
		assertEquals("chr1", intersectionList.get(2).getChr());
		assertEquals("chr2", intersectionList.get(3).getChr());
		assertEquals("chr2", intersectionList.get(4).getChr());
		assertEquals("chr2", intersectionList.get(5).getChr());
		assertEquals(5, intersectionList.get(0).getStart());
		assertEquals(10, intersectionList.get(1).getStart());
		assertEquals(15, intersectionList.get(2).getStart());
		assertEquals(5, intersectionList.get(3).getStart());
		assertEquals(10, intersectionList.get(4).getStart());
		assertEquals(15, intersectionList.get(5).getStart());
		assertEquals(9, intersectionList.get(0).getEnd());
		assertEquals(14, intersectionList.get(1).getEnd());
		assertEquals(20, intersectionList.get(2).getEnd());
		assertEquals(9, intersectionList.get(3).getEnd());
		assertEquals(14, intersectionList.get(4).getEnd());
		assertEquals(20, intersectionList.get(5).getEnd());
		assertEquals("1-1", intersectionList.get(0).getName());
		assertEquals("1-2", intersectionList.get(1).getName());
		assertEquals("1-3", intersectionList.get(2).getName());
		assertEquals("1-1", intersectionList.get(3).getName());
		assertEquals("1-2", intersectionList.get(4).getName());
		assertEquals("1-3", intersectionList.get(5).getName());

		assertEquals(3, a1.getIntersectedRegions().size());
		assertTrue(a1.getIntersectedRegions().contains(b1));
		assertTrue(a1.getIntersectedRegions().contains(b2));
		assertTrue(a1.getIntersectedRegions().contains(b3));
		assertTrue(b1.getIntersectedRegions().contains(a1));
		assertTrue(b2.getIntersectedRegions().contains(a1));
		assertTrue(b3.getIntersectedRegions().contains(a1));

		assertEquals(3, a1_2.getIntersectedRegions().size());
		assertTrue(a1_2.getIntersectedRegions().contains(b1_2));
		assertTrue(a1_2.getIntersectedRegions().contains(b2_2));
		assertTrue(a1_2.getIntersectedRegions().contains(b3_2));
		assertTrue(b1_2.getIntersectedRegions().contains(a1_2));
		assertTrue(b2_2.getIntersectedRegions().contains(a1_2));
		assertTrue(b3_2.getIntersectedRegions().contains(a1_2));
	}
}