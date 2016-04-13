package edu.cwru.cbc.ASM.detect;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

import static org.testng.AssertJUnit.assertEquals;

public class FDRControlTest {
	List<Double> pValueList_normal;
	List<Double> pValueList_all;
	List<Double> pValueList_none;

	@BeforeMethod
	public void setUp() throws Exception {
		pValueList_normal = Arrays.asList(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459,
				0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.000);
		pValueList_all = Arrays.asList(0.001, 0.0011, 0.00111, 0.001111);
		pValueList_none = Arrays.asList(0.2, 0.3, 0.4, 0.5);
	}

	@Test
	public void testGetBHYFDRCutoff() throws Exception {
		assertEquals(0.0019, FDRControl.getBHYFDRCutoff(pValueList_normal, 0.05), 0.0000001);
		assertEquals(0.001111, FDRControl.getBHYFDRCutoff(pValueList_all, 0.05), 0.0000001);
		assertEquals(-1.0, FDRControl.getBHYFDRCutoff(pValueList_none, 0.05), 0.0000001);
	}

	@Test
	public void testGetBHFDRCutoff() throws Exception {
		assertEquals(0.0095, FDRControl.getBHFDRCutoff(pValueList_normal, 0.05), 0.0000001);
		assertEquals(0.001111, FDRControl.getBHFDRCutoff(pValueList_all, 0.05), 0.0000001);
		assertEquals(-1.0, FDRControl.getBHFDRCutoff(pValueList_none, 0.05), 0.0000001);
	}
}