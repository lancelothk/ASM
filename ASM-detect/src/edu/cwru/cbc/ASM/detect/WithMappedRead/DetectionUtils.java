package edu.cwru.cbc.ASM.detect.WithMappedRead;

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
     *
     * @param pvalueList p value list
     * @param fdr        FDR
     * @return p value cutoff.
     */
    public static double getFDRCutoff(List<Double> pvalueList, double fdr) {
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
        throw new RuntimeException("no P value <= k/m*fdr!");
    }

    public static void main(String[] args) {
        //        Double[] ps = {0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.000};
        //        System.out.println(getFDRCutoff(Arrays.asList(ps), 0.05));
    }
}
