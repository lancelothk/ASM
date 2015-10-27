package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by kehu on 10/22/15.
 * Program to check if SNP is consistent with Grouping.
 */
public class SNPCheckPgm {
	private static final Splitter tabSplitter = Splitter.on("\t");

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input grouped read file").build());
		options.addOption(Option.builder("p").hasArg().desc("SNP position").build());
		options.addOption(Option.builder("a").hasArg().desc("allele pair, e.g. A-G").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String groupedReadFile = cmd.getOptionValue("i");
		int snpPosition = Integer.parseInt(cmd.getOptionValue("p"));
		String allelePair = cmd.getOptionValue("a");
		char allele1 = allelePair.charAt(0);
		char allele2 = allelePair.charAt(2);

		String[] items = groupedReadFile.replace(".mappedreads.groups.aligned", "").split("-");
		if (items.length != 3) {
			throw new RuntimeException("invalid input file name format!\t" + groupedReadFile);
		}
		int startPos = Integer.parseInt(items[1]);

		Files.readLines(new File(groupedReadFile), Charsets.UTF_8,
				new LineProcessor() {
					int[][][] observation = new int[2][2][5];
					int group = 0;
					String ref;
					int snpIndex = snpPosition - startPos;

					@Override
					public boolean processLine(@Nonnull String line) throws IOException {
						List<String> itemList = tabSplitter.splitToList(line);
						if (line.startsWith("ref:")) {
							ref = itemList.get(1);
						} else if (line.equals("")) {
							group++;
						} else {
							if (group == 2) {
								throw new RuntimeException("more than 2 groups!");
							}
							if (itemList.get(1).equals("+")) {
								checkSNPPosition(itemList.get(4), 0);
							} else if (itemList.get(1).equals("-")) {
								checkSNPPosition(itemList.get(4), 1);
							} else {
								throw new RuntimeException("invalid strand!");
							}
						}
						return true;
					}

					private void checkSNPPosition(String sequence, int strand) {
						if (snpIndex > 0 && snpIndex < sequence.length()) {
							char c = sequence.charAt(snpIndex);
							switch (c) {
								case 'A':
									observation[group][strand][0]++;
									break;
								case 'C':
									observation[group][strand][1]++;
									break;
								case 'G':
									observation[group][strand][2]++;
									break;
								case 'T':
									observation[group][strand][3]++;
									break;
								case '.':
									break;
								default:
									// 'N'
									observation[group][strand][4]++;
							}
						}
					}


					private Set<Character> plusStrandExpectedAllele(char allele) {
						Set<Character> expectedSNPSet = new HashSet<>();
						switch (allele) {
							case 'A':
								expectedSNPSet.add('A');
								break;
							case 'C':
								if (snpIndex < ref.length() - 1 && ref.charAt(snpIndex + 1) == 'G') {
									// is in CpG
									expectedSNPSet.add('C');
									expectedSNPSet.add('T');
								} else {
									expectedSNPSet.add('T');
								}
								break;
							case 'G':
								expectedSNPSet.add('G');
								break;
							case 'T':
								expectedSNPSet.add('T');
								break;
							default:
								throw new RuntimeException("invalid allele!");
						}
						return expectedSNPSet;
					}

					private Set<Character> minusStrandExpectedAllele(char allele) {
						Set<Character> expectedSNPSet = new HashSet<>();
						switch (allele) {
							case 'A':
								expectedSNPSet.add('A');
								break;
							case 'C':
								expectedSNPSet.add('C');
								break;
							case 'G':
								if (snpIndex > 1 && ref.charAt(snpIndex - 1) == 'C') {
									// is in CpG
									expectedSNPSet.add('G');
									expectedSNPSet.add('A');
								} else {
									expectedSNPSet.add('A');
								}
								break;
							case 'T':
								expectedSNPSet.add('T');
								break;
							default:
								throw new RuntimeException("invalid allele!");
						}
						return expectedSNPSet;
					}

					@Override
					public Boolean getResult() {
						// first dimension is group
						// second dimension is strand
						Map[][] observedAlleles = new Map[2][2];
						Set[][] expectedAlleles = new Set[2][2];
						for (int i = 0; i < 2; i++) {
							System.out.println("group" + (i + 1));
							System.out.println("  A\tC\tG\tT\tN");
							System.out.println(String.format("+:%d\t%d\t%d\t%d\t%d", observation[i][0][0],
									observation[i][0][1], observation[i][0][2], observation[i][0][3],
									observation[i][0][4]));
							System.out.println(String.format("-:%d\t%d\t%d\t%d\t%d", observation[i][1][0],
									observation[i][1][1], observation[i][1][2], observation[i][1][3],
									observation[i][1][4]));

							observedAlleles[i][0] = observedAlleles(observation[i][0]);
							observedAlleles[i][1] = observedAlleles(observation[i][1]);
							System.out.printf("Observed alleles:%s\t%s\n", observedAlleles[i][0],
									observedAlleles[i][1]);
						}
						expectedAlleles[0][0] = plusStrandExpectedAllele(allele1);
						expectedAlleles[0][1] = minusStrandExpectedAllele(allele1);
						expectedAlleles[1][0] = plusStrandExpectedAllele(allele2);
						expectedAlleles[1][1] = minusStrandExpectedAllele(allele2);
						System.out.printf("expected allele pair for allele1:%s\t%s\n", expectedAlleles[0][0]
								, expectedAlleles[0][1]);
						System.out.printf("expected allele pair for allele2:%s\t%s\n", expectedAlleles[1][0],
								expectedAlleles[1][1]);

						boolean isConsistent = ((isConsistent(observedAlleles[0][0],
								expectedAlleles[0][0]) && isConsistent(
								observedAlleles[0][1], expectedAlleles[0][1]))  // observed allele 1 ~ expected allele 1
								&& (isConsistent(observedAlleles[1][0], expectedAlleles[1][0]) && isConsistent(
								observedAlleles[1][1], expectedAlleles[1][1]))) // observed allele 2 ~ expected allele 2
								||
								((isConsistent(observedAlleles[0][0], expectedAlleles[1][0]) && isConsistent(
										observedAlleles[0][1],
										expectedAlleles[1][1])) // observed allele 1 ~ expected allele 2
										&& (isConsistent(observedAlleles[1][0], expectedAlleles[0][0]) && isConsistent(
										observedAlleles[1][1],
										expectedAlleles[0][1]))); // observed allele 2 ~ expected allele 1

						if (isConsistent) {
							System.out.println("SNP is consistent with grouping!");
						} else {
							System.out.println("SNP is inconsistent with grouping!");
						}

						boolean isMajorityConsistent = ((isMajorityConsistent(observedAlleles[0][0],
								expectedAlleles[0][0]) && isMajorityConsistent(
								observedAlleles[0][1], expectedAlleles[0][1]))  // observed allele 1 ~ expected allele 1
								&& (isMajorityConsistent(observedAlleles[1][0],
								expectedAlleles[1][0]) && isMajorityConsistent(
								observedAlleles[1][1], expectedAlleles[1][1]))) // observed allele 2 ~ expected allele 2
								||
								((isMajorityConsistent(observedAlleles[0][0],
										expectedAlleles[1][0]) && isMajorityConsistent(
										observedAlleles[0][1],
										expectedAlleles[1][1])) // observed allele 1 ~ expected allele 2
										&& (isMajorityConsistent(observedAlleles[1][0],
										expectedAlleles[0][0]) && isMajorityConsistent(
										observedAlleles[1][1],
										expectedAlleles[0][1]))); // observed allele 2 ~ expected allele 1

						if (isMajorityConsistent) {
							System.out.println(
									"Majority SNP is consistent with grouping!" + ((allele1 > allele2) ? allele1 + "" + allele2 : allele2 + "" + allele1));
						} else {
							System.out.println(
									"Majority SNP is inconsistent with grouping!" + ((allele1 > allele2) ? allele1 + "" + allele2 : allele2 + "" + allele1));
						}
						return true;
					}

					private boolean isConsistent(Map observed, Set expected) {
						int inconsistentCount = 0;
						@SuppressWarnings("unchecked")
						Set<Character> expectedSNP = expected;
						@SuppressWarnings("unchecked")
						Map<Character, Integer> observedSNP = observed;
						for (Map.Entry<Character, Integer> entry : observedSNP.entrySet()) {
							if (!expectedSNP.contains(entry.getKey())) {
								inconsistentCount += entry.getValue();
							}
						}
						return inconsistentCount == 0;
					}

					private boolean isMajorityConsistent(Map observed, Set expected) {
						int consistentCount = 0, inconsistentCount = 0;
						@SuppressWarnings("unchecked")
						Set<Character> expectedSNP = expected;
						@SuppressWarnings("unchecked")
						Map<Character, Integer> observedSNP = observed;
						for (Map.Entry<Character, Integer> entry : observedSNP.entrySet()) {
							if (expectedSNP.contains(entry.getKey())) {
								consistentCount += entry.getValue();
							} else {
								inconsistentCount += entry.getValue();
							}
						}
//						return inconsistentCount == 0;
						return consistentCount >= inconsistentCount;
					}

					private Map<Character, Integer> observedAlleles(int[] obs) {
						Map<Character, Integer> observedSNPMap = new HashMap<>();
						for (int i = 0; i < obs.length; i++) {
							if (obs[i] != 0) {
								switch (i) {
									case 0:
										observedSNPMap.put('A', obs[i]);
										break;
									case 1:
										observedSNPMap.put('C', obs[i]);
										break;
									case 2:
										observedSNPMap.put('G', obs[i]);
										break;
									case 3:
										observedSNPMap.put('T', obs[i]);
										break;
									default:
								}
							}
						}
						return observedSNPMap;
					}

				});
	}
}