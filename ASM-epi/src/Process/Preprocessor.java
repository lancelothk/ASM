package Process;

import DataType.*;
import com.google.common.base.Charsets;
import com.google.common.io.Files;

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
public class Preprocessor {
	public static void main(String[] args) throws IOException {
		splitEpigenome(args[0], args[1], args[2]);
//		splitEpigenome("/home/lancelothk/chr22.fa", "/home/lancelothk/h1_r1chr22", "/home/lancelothk/test");
	}

	public static void splitEpigenome(String referenceGenomeFileName, String mappedReadFileName,
									  String outputPath) throws IOException {
		long start = System.currentTimeMillis();
		RefChr refChr = readReferenceGenome(referenceGenomeFileName);
		Map<Integer, CpGSite> refMap = extractCpGSite(refChr.getRefString());
		System.out.println("load refMap complete");
		List<MappedRead> mappedReadList = readMappedReads(refMap, mappedReadFileName);
		System.out.println("load mappedReadList complete");

		// split regions by continuous CpG coverage
		List<CpGSite> refCpGSites = new ArrayList<>(refMap.values());
		Collections.sort(refCpGSites, (CpGSite c1, CpGSite c2) -> c1.getPos() - c2.getPos());
		List<List<CpGSite>> cpgSiteIntervalList = new ArrayList<>();
		boolean cont = true;
		List<CpGSite> cpGSiteList = new ArrayList<>();
		for (CpGSite refCpGSite : refCpGSites) {
			if (refCpGSite.getCoverage() == 0) {
				cont = false;
			} else {
				if (!cont) {
					cpgSiteIntervalList.add(cpGSiteList);
					cpGSiteList = new ArrayList<>();
				}
				cpGSiteList.add(refCpGSite);
				cont = true;
			}
		}
		System.out.println("Interval count:\t" + cpgSiteIntervalList.size());
		BufferedWriter intervalSummaryWriter = new BufferedWriter(new FileWriter(
				String.format("%s/%s-intervalSummary", outputPath, refChr.getChr())));
		intervalSummaryWriter.write("chr\tstart\tend\treadCount\tCpGCount\n");
		cpgSiteIntervalList.forEach((list) -> {
			Set<MappedRead> mappedReadSet = new HashSet<MappedRead>();
			list.forEach((cpGSite) -> {
				cpGSite.getCpGList().forEach((CpG cpg) -> {
					mappedReadSet.add(cpg.getMappedRead());
				});
			});
			int startPos = mappedReadSet.stream().min((r1, r2)->r1.getStart() - r2.getStart()).get().getStart();
			int endPos = mappedReadSet.stream().max((r1, r2) -> r1.getStart() - r2.getStart()).get().getEnd();
			try {
				intervalSummaryWriter.write(String.format("%s\t%d\t%d\t%d\t%d\n", refChr.getChr(), startPos,
														  endPos, mappedReadSet.size(),
														  list.size()));
				BufferedWriter mappedReadWriter = new BufferedWriter(new FileWriter(
						String.format("%s/%s-%d-%d", outputPath, refChr.getChr(), startPos,
									  endPos)));
				mappedReadWriter.write(String.format("ref:\t%s\n", refChr.getRefString().substring(startPos, endPos+1)));
				for (MappedRead mappedRead : mappedReadSet) {
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
		});
		intervalSummaryWriter.close();
		long end = (System.currentTimeMillis() - start);
		System.out.println(end + "ms");
	}

	private static List<MappedRead> readMappedReads(Map<Integer, CpGSite> refMap,
													String mappedReadFileName) throws IOException {
		return Files.readLines(new File(mappedReadFileName), Charsets.UTF_8, new MappedReadLineProcessor(refMap));
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
