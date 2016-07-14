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
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String groupedReadFile = cmd.getOptionValue("i");

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
			List<RefCpG> refCpGList = extractCpGSite(result.getLeft(), minStart);
			HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
			for (RefCpG refCpG : refCpGList) {
				refMap.put(refCpG.getPos(), refCpG);
			}
			for (MappedRead mappedRead : group) {
				mappedRead.generateCpGsInRead(refMap);
			}
			System.out.printf("group %d:\t", groupIndex);
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
			groupIndex++;
		}
	}
}
