package edu.cwru.cbc.ASM.detect.WithMappedRead.DataType;

import com.google.common.collect.ImmutableList;
import com.sun.deploy.util.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 3/27/15.
 * <p>
 * Class for format interval detection summary
 */
public class IntervalDetectionSummary {
    private static int elementCount;
    private static String headLine;
    private static String formatLine;
    private List<Object> argumentList;
    private String summaryString;
    private double regionP;


    public IntervalDetectionSummary() {
        argumentList = new ArrayList<>();
    }

    // TODO implement argument builder here.

    public static void initializeFormat(ImmutableList<Pair<String, String>> elementPairList) {
        elementCount = elementPairList.size();
        ArrayList<String> headLineList = new ArrayList<>();
        ArrayList<String> formatLineList = new ArrayList<>();
        elementPairList.stream()
                .forEach(e -> {
                    headLineList.add(e.getLeft());
                    formatLineList.add(e.getRight());
                });
        headLine = buildLine(headLineList);
        formatLine = buildLine(formatLineList);
    }

    private static String buildLine(ArrayList<String> formatLine) {
        StringBuilder formatLineBuilder = new StringBuilder(StringUtils.join(formatLine, "\t"));
        formatLineBuilder.replace(formatLineBuilder.length() - 1, formatLineBuilder.length(), "\n");
        return formatLineBuilder.toString();
    }

    public static String getHeadLine() {
        if (headLine == null) {
            throw new RuntimeException("uninitialized head line!");
        }
        return headLine;
    }

    public void appendArgument(Object argument) {
        argumentList.add(argument);
    }

    public double getRegionP() {
        return regionP;
    }

    public void setRegionP(double regionP) {
        this.regionP = regionP;
    }

    public String getSummaryString() {
        if (elementCount != argumentList.size()) {
            throw new RuntimeException("unmatched number of elements in arguments list!");
        }
        this.summaryString = String.format(formatLine, argumentList);
        return summaryString;
    }
}
