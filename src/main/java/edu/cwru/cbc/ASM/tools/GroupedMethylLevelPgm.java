package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.io.GroupedReadsLineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.Pair;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by kehu on 7/13/16.
 */
public class GroupedMethylLevelPgm {
	private static final int ALIGN_COL_SIZE = 12;

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().required().desc("input grouped read file").build());
		options.addOption(Option.builder("p").hasArg().desc("SNP position").build());
		options.addOption(Option.builder("a").hasArg().desc("allele pair, e.g. A-G").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String groupedReadFile = cmd.getOptionValue("i");
		int snpPosition = -1;
		char allele1 = ' ';
		char allele2 = ' ';
		if (cmd.hasOption("p") && cmd.hasOption("a")) {
			snpPosition = Integer.parseInt(cmd.getOptionValue("p"));
			String allelePair = cmd.getOptionValue("a");
			allele1 = allelePair.charAt(0);
			allele2 = allelePair.charAt(2);
		}

		Pair<String, List<List<MappedRead>>> result = Files.readLines(new File(groupedReadFile),
				Charsets.UTF_8, new GroupedReadsLineProcessor());
		int minStart = result.getRight()
				.stream()
				.flatMapToInt(l -> l.stream().mapToInt(MappedRead::getStart))
				.min()
				.getAsInt();

		int groupIndex = 1;
		for (List<MappedRead> group : result.getRight()
				.stream()
				.sorted((l1, l2) -> Integer.compare(l2.size(), l1.size()))
				.collect(
						Collectors.toList())) {
			System.out.printf("group %d:\t", groupIndex);
			printGroupStat(result.getLeft(), minStart, group);
			if (snpPosition != -1 && allele1 != ' ' && allele2 != ' ') {
				List<MappedRead> group1 = new ArrayList<>();
				List<MappedRead> group2 = new ArrayList<>();
				for (MappedRead mappedRead : group) {
					if (mappedRead.getSequence()
							.length() <= snpPosition - mappedRead.getStart()) { // sequence didn't cover snp
						continue;
					}
					if (mappedRead.getStrand() == '+') {
						if (mappedRead.getSequence().charAt(snpPosition - mappedRead.getStart()) == allele1) {
							group1.add(mappedRead);
						} else if (mappedRead.getSequence()
								.charAt(snpPosition - mappedRead.getStart()) == allele2) {
							group2.add(mappedRead);
						}
					} else if (mappedRead.getStrand() == '-') {
						if (mappedRead.getSequence()
								.charAt(snpPosition - mappedRead.getStart()) == getComplementaryNucleotide(
								allele1)) {
							group1.add(mappedRead);
						} else if (mappedRead.getSequence()
								.charAt(snpPosition - mappedRead.getStart()) == getComplementaryNucleotide(
								allele2)) {
							group2.add(mappedRead);
						}
					} else {
						throw new RuntimeException("unknown strand type!");
					}

				}
				System.out.print("allele 1:\t");
				printGroupStat(result.getLeft(), minStart, group1);
				System.out.print("allele 2:\t");
				printGroupStat(result.getLeft(), minStart, group2);
			}
			groupIndex++;
		}
	}

	private static char getComplementaryNucleotide(char plus) {
		switch (plus) {
			case 'A':
				return 'T';
			case 'C':
				return 'G';
			case 'G':
				return 'C';
			case 'T':
				return 'A';
			default:
				throw new RuntimeException("unknown nucleotide character!");
		}
	}

	private static void printGroupStat(String reference, int minStart, List<MappedRead> group) {
		List<RefCpG> refCpGList = extractCpGSite(reference, minStart);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		for (MappedRead mappedRead : group) {
			mappedRead.generateCpGsInRead(refMap);
		}
		for (RefCpG refCpG : refCpGList) {
			if (refCpG.getCoveredCount() != 0) {
				System.out.print(Strings.padEnd(
						String.format("%.2f(%d)", refCpG.getMethylLevel(), refCpG.getCoveredCount()),
						ALIGN_COL_SIZE, ' '));
			} else {
				// fill the gap
				System.out.print(Strings.repeat(" ", ALIGN_COL_SIZE));
			}
		}
		System.out.println();
	}
}
