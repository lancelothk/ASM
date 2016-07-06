package edu.cwru.cbc.ASM.tools;

import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;

/**
 * Created by kehu on 7/6/16.
 */
public class StrandBiasPgmTest {
	@Test
	public void testCalculateStrandBias() throws Exception {
		assertEquals(StrandBiasPgm.calculateStrandBias(11, 2, 20, 0), 2.54, 0.01);
		assertEquals(StrandBiasPgm.calculateStrandBias(16, 2, 10, 0), 1.56, 0.01);
		assertEquals(StrandBiasPgm.calculateStrandBias(8, 2, 16, 0), 2.6, 0.01);

		assertEquals(StrandBiasPgm.calculateGATKStrandBias(11, 2, 20, 0), 0.16, 0.01);
		assertEquals(StrandBiasPgm.calculateGATKStrandBias(16, 2, 10, 0), 0.12, 0.01);
		assertEquals(StrandBiasPgm.calculateGATKStrandBias(8, 2, 16, 0), 0.21, 0.01);

		assertEquals(StrandBiasPgm.calculateFisherStrandBias(11, 2, 20, 0), 0.85, 0.01);
		assertEquals(StrandBiasPgm.calculateFisherStrandBias(16, 2, 10, 0), 0.48, 0.01);
		assertEquals(StrandBiasPgm.calculateFisherStrandBias(8, 2, 16, 0), 0.86, 0.01);
	}
}