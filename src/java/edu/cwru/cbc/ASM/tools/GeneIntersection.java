package edu.cwru.cbc.ASM.tools;


import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by kehu on 6/30/15.
 */
public class GeneIntersection {
	public static void main(String[] args) throws IOException {
		String fileB = "/home/lancelothk/experiments/ASM/analysis/imprintGene/WAMIDEX/WAMIDEX.gene";
		intersectGenes("/home/lancelothk/experiments/ASM/analysis/original_gene/h1_refSeq.gene",
				fileB);
		intersectGenes("/home/lancelothk/experiments/ASM/analysis/original_gene/h9_refSeq.gene",
				fileB);
		intersectGenes("/home/lancelothk/experiments/ASM/analysis/original_gene/i90_refSeq.gene",
				fileB);
		intersectGenes("/home/lancelothk/experiments/ASM/analysis/original_gene/i90_ipsc_refSeq.gene",
				fileB);

//		notIntersectGenes("/home/lancelothk/experiments/ASM/analysis/imprintGene/WAMIDEX/WAMIDEX.gene","/home/lancelothk/WAMIDEX");
	}

	private static void notIntersectGenes(String fileA, String fileB) throws IOException {
		List<String> listA = Files.readLines(
				new File(fileA),
				Charsets.UTF_8, new LineReader());
		List<String> listB = Files.readLines(
				new File(fileB), Charsets.UTF_8,
				new LineReader());
		HashSet<String> resultSet = new HashSet<>();
		for (String sA : listA) {
			boolean found = false;
			for (String sB : listB) {
				if (sA.equals(sB)) {
					found = true;
					break;
				}
			}
			if (!found) {
				resultSet.add(sA);
			}
		}
		resultSet.forEach(System.out::println);
		System.out.println(fileA + "\t" + resultSet.size());
	}

	private static void intersectGenes(String fileA, String fileB) throws IOException {
		List<String> listA = Files.readLines(
				new File(fileA),
				Charsets.UTF_8, new LineReader());
		List<String> listB = Files.readLines(
				new File(fileB), Charsets.UTF_8,
				new LineReader());
		HashSet<String> resultSet = new HashSet<>();
		for (String sA : listA) {
			for (String sB : listB) {
				if (sA.equals(sB)) {
					resultSet.add(sB);
				}
			}
		}
		resultSet.forEach(System.out::println);
		System.out.println(fileA + "\t" + resultSet.size());
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
