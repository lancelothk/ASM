package edu.cwru.cbc.ASM.simulation;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.commons.CommonsUtils;
import edu.cwru.cbc.ASM.commons.CpG.CpG;
import edu.cwru.cbc.ASM.commons.CpG.RefChr;
import edu.cwru.cbc.ASM.commons.CpG.RefCpG;
import edu.cwru.cbc.ASM.commons.GenomicInterval;
import edu.cwru.cbc.ASM.commons.MethylStatus;
import edu.cwru.cbc.ASM.commons.Read.MappedRead;
import edu.cwru.cbc.ASM.commons.Read.MappedReadLineProcessor;
import edu.cwru.cbc.ASM.commons.bed.BedUtils;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/12/15.
 * generate semi-simulated ASM data
 */
public class Simulation {
	public static void main(String[] args) throws IOException, ParseException {
		long start = System.currentTimeMillis();

		Options options = new Options();
		options.addOption("r", true, "Reference File");
		options.addOption("m", true, "MappedRead File");
		options.addOption("t", true, "Target region File");
		options.addOption("o", true, "Output file");
		options.addOption("a", true, "Alpha");

		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse(options, args);

		String referenceGenomeFileName = cmd.getOptionValue("r");
		String readsFileName = cmd.getOptionValue("m");
		String targetRegionFileName = cmd.getOptionValue("t");
		String outputFileName = cmd.getOptionValue("o");
		double alpha = Double.valueOf(cmd.getOptionValue("a"));

		executeSimulation(referenceGenomeFileName, readsFileName, targetRegionFileName, outputFileName, alpha);
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");
	}

	/**
	 * execute data simulation
	 *
	 * @param alpha read will keep same pattern with chosen allele with probability alpha.
	 * @throws IOException
	 */
	public static void executeSimulation(String referenceGenomeFileName, String readsFileName,
			String targetRegionFileName, String outputFileName, double alpha) throws IOException {
		File outputFile = new File(outputFileName);
		if (!outputFile.getParentFile().exists()) {
			if (!outputFile.getParentFile().mkdirs()) {
				throw new RuntimeException("failed to mkdir:" + outputFile.getParentFile());
			}
		}

		System.out.println("simulation start");
		// read reference and refCpGs
		RefChr refChr = CommonsUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = CommonsUtils.extractCpGSite(refChr.getRefString(), 0);

		// read target regions
		List<GenomicInterval> targetRegionsMap = BedUtils.readSingleChromBedRegions(targetRegionFileName);
		Collections.sort(targetRegionsMap);

		// generate non-ASM regions
		List<GenomicInterval> nonASMRegions = generateNonASMRegions(refChr, targetRegionsMap);

		// attach refCpG to regions and generate random allele pattern.
		attachRefCpGToRegions(refCpGList, targetRegionsMap, nonASMRegions);

		// read input sequences
		List<MappedRead> mappedReadList = Files.asCharSource(new File(readsFileName), Charsets.UTF_8)
				.readLines(new MappedReadLineProcessor(refCpGList, 0));

		System.out.println("load reads finished");

		// assign methyl status
		assignMethylStatusForNonASMRegion(nonASMRegions);
		assignMethylStatusForASMRegion(targetRegionsMap);

		// add randomness
		addRandomnessToReads(nonASMRegions, alpha);
		addRandomnessToReads(targetRegionsMap, alpha);

		// write out result
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		for (MappedRead mappedRead : mappedReadList) {
			writer.write(mappedRead.toSimulationString() + "\n");
		}
		writer.close();

		BufferedWriter asmWriter = new BufferedWriter(new FileWriter(outputFileName + ".asmPattern"));
		BufferedWriter nonasmWriter = new BufferedWriter(new FileWriter(outputFileName + ".nonasmPattern"));
		for (GenomicInterval targetRegion : targetRegionsMap) {
			asmWriter.write("ASM\t" + targetRegion.toPatternString() + "\n");
		}
		for (GenomicInterval nonASMRegion : nonASMRegions) {
			nonasmWriter.write("nonASM\t" + nonASMRegion.toPatternString() + "\n");

		}
		asmWriter.close();
		nonasmWriter.close();
		System.out.println("simulation finished");
	}

	/**
	 * generate Non ASM regions which covers the region not covered by ASM regions.
	 *
	 * @return list of non ASM regions. In position increasing order.
	 */
	private static List<GenomicInterval> generateNonASMRegions(RefChr refChr, List<GenomicInterval> targetRegionList) {
		List<GenomicInterval> nonASMRegions = new ArrayList<>();
		int lastEnd = -1;
		for (int i = 0; i <= targetRegionList.size(); i++) {
			if (i == 0) {
				nonASMRegions.add(
						new GenomicInterval(refChr.getChr(), refChr.getStart(), targetRegionList.get(i).getStart() - 1,
								String.valueOf(nonASMRegions.size() + 1)));
				lastEnd = targetRegionList.get(i).getEnd();
			} else if (i == targetRegionList.size()) {
				nonASMRegions.add(new GenomicInterval(refChr.getChr(), lastEnd + 1, refChr.getEnd(),
						String.valueOf(nonASMRegions.size() + 1)));
			} else {
				nonASMRegions.add(
						new GenomicInterval(refChr.getChr(), lastEnd + 1, targetRegionList.get(i).getStart() - 1,
								String.valueOf(nonASMRegions.size() + 1)));
				lastEnd = targetRegionList.get(i).getEnd();
			}
		}
		return nonASMRegions;
	}

	private static void attachRefCpGToRegions(List<RefCpG> refCpGList, List<GenomicInterval> targetRegionList,
			List<GenomicInterval> nonASMRegions) {
		List<GenomicInterval> allRegions = new ArrayList<>();
		allRegions.addAll(nonASMRegions);
		allRegions.addAll(targetRegionList);
		for (RefCpG refCpG : refCpGList) {
			for (GenomicInterval region : allRegions) {
				if (refCpG.getPos() >= region.getStart() && refCpG.getPos() <= region.getEnd()) {
					region.addRefCpG(refCpG);
					break;
				}
			}
		}
		allRegions.forEach(GenomicInterval::generateRandomAllelePattern);
	}

	private static void assignMethylStatusForNonASMRegion(List<GenomicInterval> nonASMRegions) {
		// each refCpG is attached to a region, each CpG in read belongs to a refCpG.
		// Non ASM region
		for (GenomicInterval nonASMRegion : nonASMRegions) {
			for (int i = 0; i < nonASMRegion.getRefCpGList().size(); i++) {
				for (CpG cpg : nonASMRegion.getRefCpGList().get(i).getCpGList()) {
					if (nonASMRegion.getRefMethylStatus(i)) {
						cpg.setMethylStatus(MethylStatus.C);
					} else {
						cpg.setMethylStatus(MethylStatus.T);
					}
				}
			}
		}
	}

	private static void assignMethylStatusForASMRegion(List<GenomicInterval> targetRegionList) {
		// each refCpG is attached to a region, each CpG in read belongs to a refCpG.
		// ASM region
		for (GenomicInterval targetRegion : targetRegionList) {
			Set<MappedRead> mappedReadSet = getReadsInRegion(targetRegion);

			// randomly assign read to each allele
			Random rand = new Random();
			for (MappedRead mappedRead : mappedReadSet) {
				if (rand.nextBoolean()) {
					// allele 1
					mappedRead.getCpgList()
							.stream()
							.filter(cpg -> targetRegion.getRefCpGList().contains(cpg.getRefCpG()))
							.forEach(cpg -> {
								if (targetRegion.getRefMethylStatus(cpg.getRefCpG().getIndex())) {
									cpg.setMethylStatus(MethylStatus.C);
								} else {
									cpg.setMethylStatus(MethylStatus.T);
								}
							});
				} else {
					// allele 2 : counterpart of allele 1
					mappedRead.getCpgList()
							.stream()
							.filter(cpg -> targetRegion.getRefCpGList().contains(cpg.getRefCpG()))
							.forEach(cpg -> {
								if (targetRegion.getRefMethylStatus(cpg.getRefCpG().getIndex())) {
									cpg.setMethylStatus(MethylStatus.T); // opposite to allele 1
								} else {
									cpg.setMethylStatus(MethylStatus.C); // opposite to allele 1
								}
							});
				}
			}
		}
	}

	private static void addRandomnessToReads(List<GenomicInterval> regions, double alpha) {
		for (GenomicInterval region : regions) {
			Set<MappedRead> readsInRegion = getReadsInRegion(region);
			for (MappedRead read : readsInRegion) {
				// add randomness. Flip methyl status
				read.getCpgList().stream().forEach(cpG -> {
					if (Math.random() < alpha) {
						// add randomness. Flip methyl status
						cpG.setMethylStatus(cpG.getMethylStatus() == MethylStatus.C ? MethylStatus.T : MethylStatus.C);
					}
				});
			}
		}
	}

	private static Set<MappedRead> getReadsInRegion(GenomicInterval region) {
		Set<MappedRead> mappedReadSet = new HashSet<>();
		for (int i = 0; i < region.getRefCpGList().size(); i++) {
			RefCpG refCpG = region.getRefCpGList().get(i);
			refCpG.assignIndex(i);
			mappedReadSet.addAll(refCpG.getCpGList().stream().map(CpG::getMappedRead).collect(Collectors.toList()));
		}
		return mappedReadSet;
	}
}
