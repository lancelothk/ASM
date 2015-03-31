package edu.cwru.cbc.ASM.detect;

import org.apache.commons.math3.special.Gamma;

/**
 * Created by kehu on 2/25/15.
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

		//	Force cell counts to be positive or zero.

		int a = Math.abs(n11);
		int b = Math.abs(n12);
		int c = Math.abs(n21);
		int d = Math.abs(n22);

		//	Compute parameters for hypergeometric
		//	distribution.
		int r = a + b;
		int s = c + d;
		int m = a + c;
		int n = b + d;
								/*	Get range of variation. */

		int lm = (0 > m - s) ? 0 : m - s;
		int um = (m < r) ? m : r;

								/*	Probability is 1 if no range of variation. */

		if ((um - lm + 2) == 0) {
			result[0] = 1.0D;
			result[1] = 1.0D;
		} else {
			//	Compute critical value of
			//	hypergeometric distribution
			//	which serves as cut-off when
			//	computing two-tailed test. *

			double cutoff = hypergeometricProbability(a, r, s, m, n);

			//	Compute Fisher's exact test values.
			//	result[0] = two-tailed value
			//	result[1] = left tail value
			//	result[2] = right tail value

			result[0] = 0.0D;
			result[1] = 0.0D;

			for (int i = lm; i <= um; i++) {
				double probability = hypergeometricProbability(i, r, s, m, n);

				if (i >= a)
					result[1] += probability;

				if (probability <= cutoff)
					result[0] += probability;
			}
			//	Make sure all values less then 1 .

			result[0] = Math.min(result[0], 1.0D);
			result[1] = Math.min(result[1], 1.0D);
		}

		return result;
	}

	protected static double hypergeometricProbability(int x, int n1d, int n2d, int nd1, int nd2) {
		int n3 = nd1 - x;
		int ndd = nd1 + nd2;

		double sum = logCombination(n1d, x) + logCombination(n2d, n3) - logCombination(ndd, nd1);

		return Math.exp(sum);
	}

	/**
	 * Compute log of number of combinations of n things taken k at a time.
	 *
	 * @param n Number of items.
	 * @param k Group size.
	 * @return Value for log of n items taken k at a time.
	 * <p>
	 * <p>
	 * log(combination(n,m))	= log(n!/m!(n-m)!)
	 * = log(n!) - log(m!) - log((n-m)!)
	 * </p>
	 */

	protected static double logCombination(int n, int k) {
		return logFactorial(n) -
				logFactorial(k) -
				logFactorial(n - k);
	}

	/**
	 * Computes the natural log of n factorial.
	 *
	 * @param n The number for which the log factorial is desired.
	 * @return The natural log of n! .
	 * <p>
	 * <p>
	 * We use the LogGamma function to compute log(n!) using the relation:
	 * </p>
	 * <p>
	 * <p>
	 * log(n!) = logGamma( n + 1 )
	 * </p>
	 */

	public static double logFactorial(int n) {
		if (n <= 1) {
			return 0.0D;
		} else {
			return Gamma.logGamma(n + 1.0D);
		}
	}

}
