package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.CPMR.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.CpG;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MappedReadLineProcessorWithSummary;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;

/**
 * Created by lancelothk on 5/26/14.
 * Main entry of the module.
 * Convert mapped reads into epigenome and split regions by continuous CpG coverage.
 */
public class CPMR {
	public static final int MIN_CONT_COVERAGE = 4;
	public static final int MIN_INTERVAL_READS = 10;
	public static final int MIN_INTERVAL_CPG = 5;
	public static int outputIntervalCount = 0;

	public static void main(String[] args) throws IOException {
		String chr = "chr20";
		String ref = "hg18_" + chr + ".fa";
		String cellLine = "i90";
		String replicate = "r1";
		String homeDirectory = System.getProperty("user.home");
		String experimentPath = String.format("%s/experiments/ASM/", homeDirectory);

		String referenceGenomeFileName = experimentPath + "/data/" + ref;
		String mappedReadFileName = String.format("%s/data/%s_%s_%s", experimentPath, cellLine, replicate, chr);
		String outputPath = String.format("%s/result_%s_%s/intervals_%s", experimentPath, cellLine, replicate, chr);

		splitEpigenome(referenceGenomeFileName, mappedReadFileName, outputPath, OutputFormat.MappredRead);
	}

	public static void splitEpigenome(String referenceGenomeFileName, String mappedReadFileName, String outputPath,
									  OutputFormat outputFormat) throws IOException {
		long start = System.currentTimeMillis();
		File outputFile = new File(outputPath);
		if (!outputFile.exists()) {
			//noinspection ResultOfMethodCallIgnored
			outputFile.mkdirs();
		}

		RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);
		System.out.println("load refMap complete");
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		refCpGList.sort((RefCpG c1, RefCpG c2) -> c1.getPos() - c2.getPos());
		System.out.println("cpg sorting complete");
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		// assign cpg id after sorted
		for (int i = 0; i < refCpGList.size(); i++) {
			refCpGList.get(i).assignIndex(i);
		}
		System.out.println("cpg id assignment complete");
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
						new MappedReadLineProcessorWithSummary(refCpGList, refChr.getRefString().length()));
		System.out.println("load mappedReadList complete");
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");

		start = System.currentTimeMillis();
		// split regions by continuous CpG coverage
		List<List<RefCpG>> cpgSiteIntervalList = new ArrayList<>();
		boolean cont = true;
		List<RefCpG> resultRefCpGList = new ArrayList<>();
		for (int i = 0; i < refCpGList.size(); i++) {
			RefCpG curr = refCpGList.get(i);
			RefCpG next = (i + 1) < refCpGList.size() ? refCpGList.get(i + 1) : null;
			// TODO check coverage and partial methylation of NEXT????
			if (curr.getCoverage() < MIN_CONT_COVERAGE && !curr.hasPartialMethyl()) {
				cont = false;
			} else {
				if (!cont) {
					cpgSiteIntervalList.add(resultRefCpGList);
					resultRefCpGList = new ArrayList<>();
				}
				resultRefCpGList.add(curr);
				cont = curr.hasCommonRead(next);
			}
		}

		writeIntervalSummary(outputPath, refChr, cpgSiteIntervalList, outputFormat);
		System.out.println("Raw Interval count:\t" + cpgSiteIntervalList.size());
		System.out.println("Output Interval count:\t" + outputIntervalCount);
		System.out.println((System.currentTimeMillis() - start) / 1000.0 + "s");
	}

	private static void writeIntervalSummary(String outputPath, RefChr refChr, List<List<RefCpG>> cpgSiteIntervalList,
											 OutputFormat outputFormat) throws IOException {
		BufferedWriter intervalSummaryWriter = new BufferedWriter(
				new FileWriter(String.format("%s/%s-intervalSummary", outputPath, refChr.getChr())));
		intervalSummaryWriter.write("chr\tstart\tend\tlength\treadCount\tCpGCount\tstartCpG\tendCpG\n");
		cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = new HashSet<>();
			list.forEach(
					(cpGSite) -> cpGSite.getCpGList().forEach((CpG cpg) -> mappedReadSet.add(cpg.getMappedRead())));
			// only pass high quality result for next step.
			if (mappedReadSet.size() >= MIN_INTERVAL_READS && list.size() >= MIN_INTERVAL_CPG) {
				outputIntervalCount++;
				int startPos = mappedReadSet.stream().min((r1, r2) -> r1.getStart() - r2.getStart()).get().getStart();
				int endPos = mappedReadSet.stream().max((r1, r2) -> r1.getStart() - r2.getStart()).get().getEnd();
				int startCpGPos = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				int endCpGPos = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
				try {
					intervalSummaryWriter.write(
							String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", refChr.getChr(), startPos, endPos,
										  endPos - startPos + 1, mappedReadSet.size(), list.size(), startCpGPos,
										  endCpGPos));
					switch (outputFormat) {
						case MappredRead:
							Utils.writeMappedReadInInterval(outputPath, refChr, startPos, endPos, mappedReadSet);
							break;
						case ExtEpiRead:
							Utils.writeExtEpireadInInterval(outputPath, refChr, startCpGPos, endCpGPos, mappedReadSet);
							break;
						default:
							throw new RuntimeException("invalid output format");
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		});
		intervalSummaryWriter.close();
	}

	private enum OutputFormat {
		MappredRead, ExtEpiRead
	}
}