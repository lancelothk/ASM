package edu.cwru.cbc.ASM.detect;

import org.apache.commons.lang3.StringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;


/**
 * Created by kehu on 12/9/14.
 * Common utils for Detection
 */
public class DetectionUtils {
	public static String readRefFromIntervalFile(File inputFile) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFile));
		String line = bufferedReader.readLine();
		String[] items = line.split("\t");
		if (items.length != 2) {
			throw new RuntimeException(
					"invalid reference line in interval read file!\t" + line + "\t" + inputFile.getName());
		} else {
			// item[1](reference string) should be in uppercase and without space
			assert StringUtils.isAllUpperCase(items[1]);
			assert !items[1].contains(" ");
			return items[1];
		}
	}

	/**
	 * Given fdr and p value list, return FDR cutoff
	 * Using Benjamini–Hochberg–Yekutieli procedure
	 *
	 * @param pvalueList p value list
	 * @param fdr        FDR alpha
	 * @return p value cutoff.
	 */
	public static double getBHYFDRCutoff(List<Double> pvalueList, double fdr) {
		pvalueList.sort(Double::compare);
		double cm = 0;
		for (double i = 1; i < pvalueList.size(); i++) {
			cm += 1 / i;
		}
		for (int i = 1; i < pvalueList.size(); i++) {
			if (pvalueList.get(i) > fdr * (i + 1) / (double) pvalueList.size() / cm) {
				System.out.println("k is " + i);
				System.out.println("n is " + pvalueList.size());
				System.out.println("cm is " + cm);
				return pvalueList.get(i - 1);
			}
		}
		throw new RuntimeException("no P value <= a*k/(m*cm)");
	}

	/**
	 * Given fdr and p value list, return FDR cutoff
	 * Using Benjamini–Hochberg procedure
	 *
	 * @param pvalueList p value list
	 * @param fdr        FDR alpha
	 * @return p value cutoff.
	 */
	public static double getBHFDRCutoff(List<Double> pvalueList, double fdr) {
		pvalueList.sort(Double::compare);
		for (int i = 1; i < pvalueList.size(); i++) {
			if (pvalueList.get(i) > fdr * (i + 1) / (double) pvalueList.size()) {
				System.out.println("k is " + i);
				System.out.println("n is " + pvalueList.size());
				return pvalueList.get(i - 1);
			}
		}
		throw new RuntimeException("no P value <= a*k/m");
	}

	public static double correctPbyBonferroni(double p, int count) {
		return p * count;
	}

	public static double correctPbySidak(double p, int count) {
		return 1 - Math.pow(1 - p, count);
	}

	//	public static List<Double> correctPbyBonferroni(List<Double> pvalueList){
	//		return pvalueList.stream().map(p->p*pvalueList.size()).collect(Collectors.toList());
	//	}
	//
	//	public static List<Double> correctPbySidak(List<Double> pvalueList){
	//		return pvalueList.stream().map(p -> 1 - Math.pow(1 - p, pvalueList.size())).collect(Collectors.toList());
	//	}

	//	public static List<Double> correctPbyBenjaminiHochberg(List<Double> pvalueList){
	//		List<Double> correctedPvalues = new ArrayList<>(pvalueList.size());
	//		for (int i = 0; i < pvalueList.size(); i++) {
	//			correctedPvalues.add(pvalueList.get(i)*pvalueList.size()/i);
	//		}
	//		return correctedPvalues;
	//	}
	//
	//	public static List<Double> correctPbyBenjaminiHochbergYekutieli(List<Double> pvalueList){
	//		List<Double> correctedPvalues = new ArrayList<>(pvalueList.size());
	//		double cm = 0;
	//		for (double i = 1; i < pvalueList.size(); i++) {
	//			cm += 1 / i;
	//		}
	//		for (int i = 0; i < pvalueList.size(); i++) {
	//			correctedPvalues.add(pvalueList.get(i)*pvalueList.size()/(i*cm));
	//		}
	//		return correctedPvalues;
	//	}
	//
	//	public static List<Double> correctPbyBonferroniHolm(List<Double> pvalueList){
	//		List<Double> correctedPvalues = new ArrayList<>(pvalueList.size());
	//		for (int i = 0; i < pvalueList.size(); i++) {
	//			correctedPvalues.add(pvalueList.get(i)*(pvalueList.size()-i+1));
	//		}
	//		return correctedPvalues;
	//	}
}
