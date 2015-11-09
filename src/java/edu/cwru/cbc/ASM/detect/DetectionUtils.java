package edu.cwru.cbc.ASM.detect;

import com.google.common.collect.Iterables;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.detect.dataType.Vertex;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Collection;
import java.util.List;

/**
 * Created by kehu on 10/20/15.
 * Utils class for ASm detection
 */

public class DetectionUtils {

	@SuppressWarnings("unused")
	public static double calcRegionP_SidakComb(List<RefCpG> twoClusterRefCpGList) {
		double minP = Double.MAX_VALUE;
		for (RefCpG refCpG : twoClusterRefCpGList) {
			if (minP > refCpG.getP_value()) {
				minP = refCpG.getP_value();
			}
		}
		return 1 - Math.pow(1 - minP, twoClusterRefCpGList.size());
	}

	@SuppressWarnings("unused")
	public static double calcRegionP_StoufferComb(List<RefCpG> twoClusterRefCpGList) {
		NormalDistribution stdNorm = new NormalDistribution(0, 1);
		double z = twoClusterRefCpGList.stream()
				.mapToDouble(refCpG -> refCpG.getP_value() == 1 ? 0 : stdNorm.inverseCumulativeProbability(
						1 - refCpG.getP_value())).sum() / Math.sqrt(twoClusterRefCpGList.size());
		return 1 - stdNorm.cumulativeProbability(z);
	}

	@SuppressWarnings("unused")
	public static double calcRegionP_WeightedStoufferComb(List<RefCpG> twoClusterRefCpGList) {
		NormalDistribution stdNorm = new NormalDistribution(0, 1);
		double z = twoClusterRefCpGList.stream()
				.mapToDouble(
						refCpG -> refCpG.getCoveredCount() * (refCpG.getP_value() == 1 ? 0 : stdNorm.inverseCumulativeProbability(
								1 - refCpG.getP_value()))).sum() / Math.sqrt(twoClusterRefCpGList.stream()
				.mapToInt(refCpG -> refCpG.getCoveredCount() * refCpG.getCoveredCount())
				.sum());
		return 1 - stdNorm.cumulativeProbability(z);
	}

	public static double calcRegionP_FisherCombSum(List<RefCpG> twoClusterRefCpGList) {
		return -2 * twoClusterRefCpGList.stream().mapToDouble(refCpG -> Math.log(refCpG.getP_value())).sum();
	}

	@SuppressWarnings("unused")
	public static double calcRegionP_FisherComb(List<RefCpG> twoClusterRefCpGList) {
		ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(2 * twoClusterRefCpGList.size());
		double combSum = -2 * twoClusterRefCpGList.stream().mapToDouble(refCpG -> Math.log(refCpG.getP_value())).sum();
		return 1 - chiSquaredDistribution.cumulativeProbability(combSum);
	}

	public static void fisherTest(List<RefCpG> twoClusterRefCpGList,
	                              Collection<Vertex> clusterResults) {
		for (RefCpG refCpG : twoClusterRefCpGList) {
			int j = 0;
			int[][] matrix = new int[2][2];
			for (Vertex vertex : clusterResults) {
				if (vertex.getRefCpGMap().containsKey(refCpG.getPos())) {
					matrix[j][0] = vertex.getRefCpGMap().get(refCpG.getPos()).getMethylCount();
					matrix[j][1] = vertex.getRefCpGMap().get(refCpG.getPos()).getNonMethylCount();
					j++;
				}
			}
			double fisher_P = FisherExactTest.fishersExactTest(matrix[0][0], matrix[0][1], matrix[1][0],
					matrix[1][1])[0];  // [0] is two tail test.
			if (fisher_P >= 0 && fisher_P <= 1) {
				refCpG.setP_value(fisher_P);
			} else {
				throw new RuntimeException("p value is not in [0,1]!");
			}

		}
	}

	public static double calcClusterIndex(List<RefCpG> twoClusterRefCpGList,
	                                      Collection<Vertex> clusterResults) {
		if (clusterResults.size() != 2) {
			return 1;
		}
		Vertex cluster1 = Iterables.get(clusterResults, 0);
		Vertex cluster2 = Iterables.get(clusterResults, 1);
		double interClusterDistance = 0, intraClusterDistance = (cluster1.getIntraClusterDistance(
				twoClusterRefCpGList) + cluster2.getIntraClusterDistance(twoClusterRefCpGList)) / 2;
		for (RefCpG refCpG : twoClusterRefCpGList) {
			interClusterDistance += Math.abs(cluster1.getRefCpGMap().get(refCpG.getPos()).getMethylLevel()
					- cluster2.getRefCpGMap().get(refCpG.getPos()).getMethylLevel());
		}
		// TODO: may consider to use Euclidean distance here
		return intraClusterDistance / interClusterDistance;
	}
}
