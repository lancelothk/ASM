package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.detect.dataType.IntervalDetectionSummaryFormatter;
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
		IntervalDetectionSummaryFormatter.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>()
						.add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("originalRegion", "%s"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#clusterRefCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("regionFScore", "%f"))
						.add(new ImmutablePair<>("clusterIndex", "%f"))
						.add(new ImmutablePair<>("#group1", "%d"))
						.add(new ImmutablePair<>("#group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.build());
	}

	@Test
	public void test_detection_FP_random() throws Exception {
		URL file = getClass().getClassLoader().getResource("chr20-25795835-25796072.mappedreads");
		assertNotNull(file);
		Detection detection = new Detection(new File(file.getFile()), ".", 5, 1000);
		String summary = detection.call();
		assertEquals("not empty summary", "", summary);
	}

	@Test
	public void test_detection_TP_H9_GANS() throws Exception {
		URL file = getClass().getClassLoader().getResource("chr20-57430213-57430368.mappedreads");
		assertNotNull(file);
		Detection detection = new Detection(new File(file.getFile()), ".", 5, 1000);
		String summary = detection.call();
		String[] itemList = summary.split("\t");
		assertEquals("MEC", "4.000000", itemList[10]);
		assertEquals("normMEC", "0.027586", itemList[11]);
		assertEquals("regionP", "0.000000e+00", itemList[12]);
		assertEquals("regionFScore", "145.027582", itemList[13]);
		assertEquals("clusterIndex", "0.003362", itemList[14]);
		assertEquals("group1 size", "13", itemList[15]);
		assertEquals("group2 size", "21", itemList[16]);
		assertEquals("group1 methyl", "0.988636", itemList[17]);
		assertEquals("group2 methyl", "0.037338", itemList[18]);
	}
}