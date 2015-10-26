package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.util.List;

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

		boolean isConsistent = Files.readLines(new File(groupedReadFile), Charsets.UTF_8,
				new LineProcessor<Boolean>() {
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


					private String plusStrandExpectedAllele(char allele) {
						switch (allele) {
							case 'A':
								return "A";
							case 'C':
								if (snpIndex < ref.length() - 1 && ref.charAt(snpIndex + 1) == 'G') {
									// is in CpG
									return "C/T";
								} else {
									return "T";
								}
							case 'G':
								return "G";
							case 'T':
								return "T";
							default:
								throw new RuntimeException("invalid allele!");
						}
					}

					private String minusStrandExpectedAllele(char allele) {
						switch (allele) {
							case 'A':
								return "A";
							case 'C':
								return "C";
							case 'G':
								if (snpIndex > 1 && ref.charAt(snpIndex - 1) == 'C') {
									// is in CpG
									return "G/A";
								} else {
									return "A";
								}
							case 'T':
								return "T";
							default:
								throw new RuntimeException("invalid allele!");
						}
					}

					@Override
					public Boolean getResult() {
						// first dimension is group
						// second dimension is strand
						String[][] observedAlleles = new String[2][2];
						String[][] expectedAlleles = new String[2][2];
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

						return ((isSubSet(observedAlleles[0][0], expectedAlleles[0][0]) && isSubSet(
								observedAlleles[0][1], expectedAlleles[0][1]))  // observed allele 1 ~ expected allele 1
								&& (isSubSet(observedAlleles[1][0], expectedAlleles[1][0]) && isSubSet(
								observedAlleles[1][1], expectedAlleles[1][1]))) // observed allele 2 ~ expected allele 2
								||
								((isSubSet(observedAlleles[0][0], expectedAlleles[1][0]) && isSubSet(
										observedAlleles[0][1],
										expectedAlleles[1][1])) // observed allele 1 ~ expected allele 2
										&& (isSubSet(observedAlleles[0][0], expectedAlleles[1][0]) && isSubSet(
										observedAlleles[1][1],
										expectedAlleles[1][1]))); // observed allele 2 ~ expected allele 1

					}

					private boolean isSubSet(String s, String sBase) {
						for (char c : s.toCharArray()) {
							if (sBase.indexOf(c) == -1) {
								return false;
							}
						}
						return true;
					}

					private String observedAlleles(int[] obs) {
						StringBuilder sb = new StringBuilder();
						for (int i = 0; i < obs.length; i++) {
							if (obs[i] != 0) {
								switch (i) {
									case 0:
										sb.append('A');
										break;
									case 1:
										sb.append('C');
										break;
									case 2:
										sb.append('G');
										break;
									case 3:
										sb.append('T');
										break;
									default:
								}
							}
						}
						return sb.toString();
					}

				});

		if (isConsistent) {
			System.out.println("SNP is consistent with grouping!");
		} else {
			System.out.println("SNP is inconsistent with grouping!");
		}
	}
}