package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.detect.dataType.IntervalDetectionSummary;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.net.URL;

import static org.testng.AssertJUnit.assertEquals;
import static org.testng.AssertJUnit.assertNotNull;

public class DetectionTest {
	@BeforeMethod
	public void setUp() throws Exception {
		IntervalDetectionSummary.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>().add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("oldRegion", "%s"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#refCpG", "%d"))
						.add(new ImmutablePair<>("#clusterCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("randPIndex", "%d"))
						.add(new ImmutablePair<>("clusterIndex", "%f"))
						.add(new ImmutablePair<>("#group1", "%d"))
						.add(new ImmutablePair<>("#group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.add(new ImmutablePair<>("FDR_label", "%s"))
						.build());
	}

	@Test
	public void test_detection() throws Exception {
		URL file = getClass().getClassLoader().getResource("chr20-25795835-25796072.mappedreads");
		assertNotNull(file);
		Detection detection = new Detection(new File(file.getFile()), 5, 6, 1000);
		IntervalDetectionSummary intervalDetectionSummary = detection.call();
		String[] itemList = intervalDetectionSummary.getSummaryString(0).split("\t");
		assertEquals("MEC", "28.000000", itemList[11]);
		assertEquals("normMEC", "0.210526", itemList[12]);
		assertEquals("regionP", "4.595401e-05", itemList[13]);
		assertEquals("dbindex", "0.047111", itemList[15]);
		assertEquals("group1 size", "21", itemList[16]);
		assertEquals("group2 size", "25", itemList[17]);
		assertEquals("group1 methyl", "0.528889", itemList[18]);
		assertEquals("group2 methyl", "0.528690", itemList[19]);
	}
}