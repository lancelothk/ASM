package edu.cwru.cbc.ASM.commons.genomicInterval;

import org.testng.annotations.Test;

import static org.testng.Assert.assertTrue;

public class GenomicIntervalBaseTest {

	@Test
	public void testCompareTo() throws Exception {
		GenomicIntervalBase interval1 = new GenomicIntervalBase("chr1", 10, 20) {
		};
		GenomicIntervalBase interval2 = new GenomicIntervalBase("chr1", 0, 20) {
		};
		GenomicIntervalBase interval3 = new GenomicIntervalBase("chr2", 0, 15) {
		};
		GenomicIntervalBase interval4 = new GenomicIntervalBase("chr2", 0, 20) {
		};

		assertTrue(interval1.compareTo(interval2) > 0); // same chr/end, diff start
		assertTrue(interval2.compareTo(interval4) < 0); // same start/end, diff chr
		assertTrue(interval3.compareTo(interval4) < 0); // same chr/start, diff end
		assertTrue(interval1.compareTo(interval3) < 0); // all diff
	}
}