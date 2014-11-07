package edu.cwru.cbc.ASM.split;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.split.DataType.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by lancelothk on 5/26/14.
 * Main entry of the module.
 * Convert mapped reads into epigenome and split regions by continuous CpG coverage.
 */
public class SplitInterval {
    public static final int MIN_CONT_COVERAGE = 1;
    public static final int MIN_INTERVAL_READS = 15;
    public static final int MIN_INTERVAL_CPG = 5;
    public static int outputCount = 0;

	public static void main(String[] args) throws IOException {
//		splitEpigenome(args[0], args[1], args[2]);

        String ref = "hg18_chr22.fa";
        String cellLine = "i90";
        String replicate = "r1";
        String experimentPath = "/home/kehu/experiments/ASM/";
        splitEpigenome(experimentPath + "/data/" + ref,
                       String.format("%s/data/%s_%s_chr22", experimentPath, cellLine, replicate),
                       experimentPath + "/test");
    }

	public static void splitEpigenome(String referenceGenomeFileName, String mappedReadFileName,
									  String outputPath) throws IOException {
		long start = System.currentTimeMillis();
		RefChr refChr = readReferenceGenome(referenceGenomeFileName);
		Map<Integer, CpGSite> refMap = extractCpGSite(refChr.getRefString());
		System.out.println("load refMap complete");
		List<MappedRead> mappedReadList = readMappedReads(refMap, mappedReadFileName, refChr.getRefString().length());
		System.out.println("load mappedReadList complete");

		// split regions by continuous CpG coverage
		List<CpGSite> refCpGSites = new ArrayList<>(refMap.values());
		Collections.sort(refCpGSites, (CpGSite c1, CpGSite c2) -> c1.getPos() - c2.getPos());
		List<List<CpGSite>> cpgSiteIntervalList = new ArrayList<>();
		boolean cont = true;
		List<CpGSite> cpGSiteList = new ArrayList<>();
		for (int i = 0; i < refCpGSites.size(); i++) {
			CpGSite curr = refCpGSites.get(i);
			CpGSite next = (i+1)<refCpGSites.size()?refCpGSites.get(i+1):null;
            if (curr.getCoverage() <= MIN_CONT_COVERAGE) {
                cont = false;
			} else {
				if (!cont) {
					cpgSiteIntervalList.add(cpGSiteList);
					cpGSiteList = new ArrayList<>();
				}
				cpGSiteList.add(curr);
				cont = hasCommonRead(curr, next);
			}
		}
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(
				String.format("%s/%s-intervalSummary", outputPath, refChr.getChr())));
        intervalSummaryWriter.write("chr\tstart\tend\tlength\treadCount\tCpGCount\tstartCpG\tendCpG\n");
        cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = new HashSet<MappedRead>();
			list.forEach((cpGSite) -> {
				cpGSite.getCpGList().forEach((CpG cpg) -> {
					mappedReadSet.add(cpg.getMappedRead());
				});
			});
            // only pass high quality result for next step.
            if (mappedReadSet.size() >= MIN_INTERVAL_READS && list.size() >= MIN_INTERVAL_CPG) {
                outputCount++;
                int startPos = mappedReadSet.stream().min((r1, r2) -> r1.getStart() - r2.getStart()).get().getStart();
				int endPos = mappedReadSet.stream().max((r1, r2) -> r1.getStart() - r2.getStart()).get().getEnd();
                int startCpG = list.stream().min((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
                int endCpG = list.stream().max((cpg1, cpg2) -> cpg1.getPos() - cpg2.getPos()).get().getPos();
                try {
                    intervalSummaryWriter.write(
                            String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", refChr.getChr(), startPos, endPos,
                                          endPos - startPos + 1, mappedReadSet.size(), list.size(), startCpG, endCpG));
                    BufferedWriter mappedReadWriter = new BufferedWriter(new FileWriter(
                            String.format("%s/%s-%d-%d", outputPath, refChr.getChr(), startPos, endPos)));
                    mappedReadWriter.write(
                            String.format("ref:\t%s\n", refChr.getRefString().substring(startPos, endPos + 1)));
                    for (MappedRead mappedRead : mappedReadSet) {
                        if (mappedRead.getCpGCount() == 0) {
                            throw new RuntimeException("interval read cover no CpG site!");
                        }
                        try {
                            mappedReadWriter.write(mappedRead.toWriteString());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                    mappedReadWriter.close();
                } catch (IOException e) {
					e.printStackTrace();
				}
			}
		});
		intervalSummaryWriter.close();
        System.out.println("Raw Interval count:\t" + cpgSiteIntervalList.size());
        System.out.println("Output Interval count:\t" + outputCount);
        long end = (System.currentTimeMillis() - start);
		System.out.println(end + "ms");
	}

	private static boolean hasCommonRead(CpGSite curr, CpGSite next) {
		for (CpG cCpg : curr.getCpGList()) {
			for (CpG nCpg : next.getCpGList()) {
				if (cCpg.getMappedRead() == nCpg.getMappedRead()){
					return true;
				}
			}
		}
		return false;
	}

	private static List<MappedRead> readMappedReads(Map<Integer, CpGSite> refMap,
													String mappedReadFileName, int refLength) throws IOException {
		return Files.readLines(new File(mappedReadFileName), Charsets.UTF_8, new MappedReadLineProcessor(refMap, refLength));
	}

	private static Map<Integer, CpGSite> extractCpGSite(String referenceGenomeString) {
		// initialize hashmap with 1/1000 of reference size
		Map<Integer, CpGSite> refMap = new HashMap<>(referenceGenomeString.length() / 1000);
		int id = 0; // 0-based id
		for (int i = 0; i < referenceGenomeString.length() - 1; i++) {
			if (referenceGenomeString.charAt(i) == 'C' && referenceGenomeString.charAt(i + 1) == 'G') {
				refMap.put(i, new CpGSite(i, id++));
			}
		}
		return refMap;
	}

	private static RefChr readReferenceGenome(String inputFileName) throws IOException {
		List<String> lines = Files.readLines(new File(inputFileName), Charsets.UTF_8);
		StringBuilder referenceBuilder = new StringBuilder();
		// parse and remove first line, which is chromosome name
		String chr = lines.get(0).replace(">", "");
		lines.remove(0);
		for (String line : lines) {
			referenceBuilder.append(line.toUpperCase());
		}
		return new RefChr(chr, referenceBuilder.toString().toUpperCase());
	}
}
