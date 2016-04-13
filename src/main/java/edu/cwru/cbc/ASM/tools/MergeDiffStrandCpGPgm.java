package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

/**
 * Created by kehu on 12/1/15.
 * Merge CpG from +/- strand into one. The input is from bismark coverage output.
 */
public class MergeDiffStrandCpGPgm {
	public static void main(String[] args) throws IOException, ParseException, ExecutionException,
			InterruptedException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input cov file").required().build());
		options.addOption(Option.builder("o").hasArg().desc("output cov file").required().build());
		options.addOption(Option.builder("a")
				.hasArg()
				.desc("alternative input cov file. Will be merged together.E.g. r1/r2")
				.build());
		options.addOption(Option.builder("r").hasArg().desc("Reference File path").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputFileName = cmd.getOptionValue("i");
		String alternativeInputFileName = cmd.getOptionValue("a", "");
		String outputFileName = cmd.getOptionValue("o");
		String referenceFilePath = cmd.getOptionValue("r");
		execution(inputFileName, alternativeInputFileName, referenceFilePath, outputFileName);
	}

	private static void execution(String inputFileName, String alternativeInputFileName, String referenceFilePath, String outputFileName) throws
			IOException {
		Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap = MethylationUtils.initializeGenomeRefCpGMap(
				IOUtils.readReferenceGenome(referenceFilePath));
		System.out.println("finished loading reference genome");
		mergeDiffStrandCpG(inputFileName, genomeRefCpGMap);
		if (!alternativeInputFileName.equals("")) {
			mergeDiffStrandCpG(alternativeInputFileName, genomeRefCpGMap);
		}
		writeSortedOutput(outputFileName, genomeRefCpGMap);
	}

	private static void writeSortedOutput(String outputFileName, Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap) throws
			IOException {
		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		List<Map.Entry<String, HashIntObjMap<RefCpG>>> entryList = new ArrayList<>(genomeRefCpGMap.entrySet());
		entryList.sort((e1, e2) -> e1.getKey().compareTo(e2.getKey()));
		for (Map.Entry<String, HashIntObjMap<RefCpG>> entry : entryList) {
			List<RefCpG> refCpGList = entry.getValue()
					.values()
					.stream()
					.sorted(RefCpG::compareTo)
					.collect(Collectors.toList());
			for (RefCpG refCpG : refCpGList) {
				bufferedWriter.write(
						String.format("%s\t%d\t%d\t%f\t%d\t%d\n", entry.getKey(), refCpG.getPos(), refCpG.getPos() + 1,
								refCpG.getMethylLevel() * 100, refCpG.getMethylCount(), refCpG.getNonMethylCount()));
			}
		}
		bufferedWriter.close();
	}

	private static void mergeDiffStrandCpG(String inputFileName, final Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap) throws
			IOException {
		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor<Object>() {
			@Override
			public boolean processLine(@Nonnull String line) throws IOException {
				List<String> items = IOUtils.tabSplitter.splitToList(line);
				// input cov is 1-based, convert to 0-based
				int pos = Integer.parseInt(items.get(1)) - 1;
				HashIntObjMap<RefCpG> chromosomeRefCpGMap = genomeRefCpGMap.get(items.get(0));
				// check if cpg come from plus strand
				RefCpG refCpG = chromosomeRefCpGMap.get(pos);
				if (refCpG != null) {
					refCpG.addMethylCount(Integer.parseInt(items.get(4)));
					refCpG.addNonMethylCount(Integer.parseInt(items.get(5)));
				} else {
					// check if cpg come from minus strand
					refCpG = chromosomeRefCpGMap.get(pos - 1);
					if (refCpG != null) {
						refCpG.addMethylCount(Integer.parseInt(items.get(4)));
						refCpG.addNonMethylCount(Integer.parseInt(items.get(5)));
					} else {
						throw new RuntimeException("Non CpG Methyl site!");
					}
				}
				return true;
			}

			@Override
			public Object getResult() {
				return null;
			}
		});

		System.out.println("finished reading input file");
	}
}
