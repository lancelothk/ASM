package edu.cwru.cbc.ASM.detect;

import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;

public class FisherExactTestTest {

	public static final double DELTA = 0.000001;

	@Test
	public void test_smallCount() throws Exception {
		assertEquals(1.0, FisherExactTest.fishersExactTest(1, 2, 3, 4)[0], DELTA);
		assertEquals(0.008284, FisherExactTest.fishersExactTest(10, 2, 2, 8)[0], DELTA);
		assertEquals(0.214286, FisherExactTest.fishersExactTest(0, 11, 3, 8)[0], DELTA);
		assertEquals(1.0, FisherExactTest.fishersExactTest(0, 0, 0, 0)[0], DELTA);
		assertEquals(0.333333, FisherExactTest.fishersExactTest(2, 0, 0, 2)[0], DELTA);
		assertEquals(0.080087, FisherExactTest.fishersExactTest(5, 1, 1, 5)[0], DELTA);
		assertEquals(0.007937, FisherExactTest.fishersExactTest(5, 0, 0, 5)[0], DELTA);
		assertEquals(0.003968, FisherExactTest.fishersExactTest(5, 0, 0, 5)[1], DELTA);

	}

	@Test
	public void test_bigCount() throws Exception {
		assertEquals(0.269604, FisherExactTest.fishersExactTest(50, 11, 40, 15)[0], DELTA);
		assertEquals(0.843795, FisherExactTest.fishersExactTest(50, 21, 40, 15)[0], DELTA);
		assertEquals(1.0, FisherExactTest.fishersExactTest(500, 500, 500, 500)[0], DELTA);
		assertEquals(0.0, FisherExactTest.fishersExactTest(1000, 0, 0, 1000)[0], DELTA);
	}
}