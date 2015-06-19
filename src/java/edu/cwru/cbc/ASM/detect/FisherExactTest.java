package edu.cwru.cbc.ASM.detect;

import org.apache.commons.math3.special.Gamma;

/**
 * Fisher's exact test.
 */
public class FisherExactTest {
	/**
	 * Calculate Fisher's exact test from the four cell counts.
	 *
	 * @return double vector with three entries.
	 * [0]	= two-sided Fisher's exact test.
	 * [1]	= right-tail Fisher's exact test.
	 */

	public static double[] fishersExactTest(int n11, int n12, int n21, int n22) {
		double result[] = new double[2];
		if (n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0) {
			throw new RuntimeException("all input should be >= 0");
		}
		int k = n11, K = n11 + n12, n = n11 + n21, N = n11 + n12 + n21 + n22;

		int mink = (n11 - n22 < 0) ? 0 : n11 - n22;
		int maxk = (n11 + n21 < K) ? n11 + n21 : K;

		/*	Probability is 1 if no range of variation. */
		if ((maxk - mink + 2) == 0) {
			result[0] = 1;
			result[1] = 1;
		} else {
			double cutoff = hyperGeometricProbability(k, K, n, N);
			for (k = mink; k <= maxk; k++) {
				double probability = hyperGeometricProbability(k, K, n, N);
				if (k >= n11) {
					result[1] += probability;
				}
				if (probability <= cutoff) {
					result[0] += probability;
				}
			}
			//	Make sure all values less then 1 .
			result[0] = Math.min(result[0], 1);
			result[1] = Math.min(result[1], 1);
		}
		return result;
	}

	private static double hyperGeometricProbability(int k, int K, int n, int N) {
		// apache common math's HyperGeometricDistribution use logBinomial to calculate probability, which is much slower than logGamma
		return Math.exp(gammaLogCombination(K, k) + gammaLogCombination(N - K, n - k) - gammaLogCombination(N, n));
	}


	private static double gammaLogCombination(int n, int k) {
		return Gamma.logGamma(n + 1.0D) - Gamma.logGamma(k + 1.0D) - Gamma.logGamma(n - k + 1.0D);
	}
}
