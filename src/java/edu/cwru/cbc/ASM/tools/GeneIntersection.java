package edu.cwru.cbc.ASM.tools;


import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 6/30/15.
 */
public class GeneIntersection {
	public static void main(String[] args) throws IOException {
		List<String> listA = Files.readLines(
				new File("/home/kehu/experiments/ASM/analysis/experimentResult/h1_annotation_refSeq_true_unique.bed"),
				Charsets.UTF_8, new LineReader());
		List<String> listB = Files.readLines(
				new File("/home/kehu/experiments/ASM/analysis/imprintGene/geneimprintcom.list"), Charsets.UTF_8,
				new LineReader());
		for (String sA : listA) {
			for (String sB : listB) {
				if (sA.equals(sB)) {
					System.out.println(sA);
				}
			}
		}
	}

	static class LineReader implements LineProcessor<List<String>> {
		private List<String> result = new ArrayList<>();

		@Override
		public boolean processLine(String line) throws IOException {
			result.add(line);
			return true;
		}

		@Override
		public List<String> getResult() {
			return result;
		}
	}
}
