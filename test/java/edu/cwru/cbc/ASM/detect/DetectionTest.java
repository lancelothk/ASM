package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.ImmutableList;
import edu.cwru.cbc.ASM.detect.dataType.IntervalDetectionSummary;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.net.URL;

import static org.testng.AssertJUnit.assertNotNull;

public class DetectionTest {
	@BeforeMethod
	public void setUp() throws Exception {
		IntervalDetectionSummary.initializeFormat(
				new ImmutableList.Builder<Pair<String, String>>().add(new ImmutablePair<>("chr", "%s"))
						.add(new ImmutablePair<>("startPos", "%d"))
						.add(new ImmutablePair<>("endPos", "%d"))
						.add(new ImmutablePair<>("length", "%d"))
						.add(new ImmutablePair<>("#edge", "%d"))
						.add(new ImmutablePair<>("#read", "%d"))
						.add(new ImmutablePair<>("#refCpG", "%d"))
						.add(new ImmutablePair<>("#clusterCpG", "%d"))
						.add(new ImmutablePair<>("#cluster", "%d"))
						.add(new ImmutablePair<>("CpGsum", "%d"))
						.add(new ImmutablePair<>("MECsum", "%f"))
						.add(new ImmutablePair<>("NormMEC", "%f"))
						.add(new ImmutablePair<>("errorProb", "%f"))
						.add(new ImmutablePair<>("regionP", "%e"))
						.add(new ImmutablePair<>("minRandP", "%e"))
						.add(new ImmutablePair<>("dbindex", "%e"))
						.add(new ImmutablePair<>("group1", "%d"))
						.add(new ImmutablePair<>("group2", "%d"))
						.add(new ImmutablePair<>("group1Methyl", "%f"))
						.add(new ImmutablePair<>("group2Methyl", "%f"))
						.add(new ImmutablePair<>("label", "%s"))
						.build());

	}

	@Test
	public void test_detection() throws Exception {
		URL file = getClass().getClassLoader().getResource("chr20-25795835-25796072.mappedreads");
		assertNotNull(file);
		Detection detection = new Detection(new File(file.getFile()), 5, 4);
		IntervalDetectionSummary intervalDetectionSummary = detection.call();
//		assertEquals(
//				"20\t25795835\t25796072\t238\t474\t46\t10\t10\t2\t133\t28.000000\t0.210526\t0.073355\t7.499676e-06\t3.772956e-02\t9.422291e-02\t21\t25\t0.528889\t0.528690\t-\n",
//				intervalDetectionSummary.getSummaryString(0));
	}
}