package edu.cwru.cbc.ASM.detect.DataType;

import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;

/**
 * Created by kehu on 3/27/15.
 * <p>
 * Class for format interval detection summary
 */
public class IntervalDetectionSummary {
	private static int elementCount;
	private static String headLine;
	private static String formatLine;
	private Object[] arguments;
	private double regionP;


	public IntervalDetectionSummary(double regionP, Object... arguments) {
		this.regionP = regionP;
		this.arguments = arguments;
	}

	public static void initializeFormat(ImmutableList<Pair<String, String>> elementPairList) {
		elementCount = elementPairList.size();
		ArrayList<String> headLineList = new ArrayList<>();
		ArrayList<String> formatLineList = new ArrayList<>();
		elementPairList.stream().forEach(e -> {
			headLineList.add(e.getLeft());
			formatLineList.add(e.getRight());
		});
		headLine = buildLine(headLineList);
		formatLine = buildLine(formatLineList);
	}

	private static String buildLine(ArrayList<String> formatLine) {
		return StringUtils.join(formatLine, "\t") + "\n";
	}

	public static String getHeadLine() {
		if (headLine == null) {
			throw new RuntimeException("uninitialized head line!");
		}
		return headLine;
	}

	public double getRegionP() {
		return regionP;
	}

	public void setRegionP(double regionP) {
		this.regionP = regionP;
	}

	public String getSummaryString(double region_threshold) {
		if (elementCount != arguments.length) {
			throw new RuntimeException(
					String.format("unmatched number of elements in arguments list!:%d-%d", elementCount,
							arguments.length + 1));
		}
		return String.format(formatLine, arguments).replace("<label>", regionP <= region_threshold ? "+" : "-");
	}

	public String getBedString() {
		return String.format("chr%s\t%d\t%d\n", arguments[0], arguments[1], arguments[2]);
	}
}
