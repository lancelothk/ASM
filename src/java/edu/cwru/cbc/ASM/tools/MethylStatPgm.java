package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.io.IOUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefChr;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.IUPACCode;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static edu.cwru.cbc.ASM.commons.methylation.MethylationUtils.extractCpGSite;

/**
 * Created by kehu on 9/16/15.
 * s
 */
public class MethylStatPgm {
	public static final int INIT_POS = 0;
	private static final Splitter tabSplitter = Splitter.on("\t");

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("c").hasArg().desc("chromosome of input").required().build());
		options.addOption(Option.builder("r").hasArg().desc("Reference File").required().build());
		options.addOption(Option.builder("m").hasArg().desc("MappedRead File").required().build());
		options.addOption(Option.builder("o").hasArg().desc("Output File").required().build());
//		options.addOption(Option.builder("p").desc("pair end mode").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String chr = cmd.getOptionValue("c");
		String referenceGenomeFileName = cmd.getOptionValue("r");
		String mappedReadFileName = cmd.getOptionValue("m");
		String outputFileName = cmd.getOptionValue("o");

		// load reference
		long start = System.currentTimeMillis();
		RefChr refChr = IOUtils.readReferenceGenome(referenceGenomeFileName);
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), INIT_POS);
		System.out.println("load refMap complete\t" + (System.currentTimeMillis() - start) / 1000.0 + "s");

		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}

		Files.readLines(new File(mappedReadFileName), Charsets.UTF_8, new LineProcessor() {
			@Override
			public boolean processLine(String line) throws IOException {
				List<String> itemList = tabSplitter.splitToList(line);
				if (!itemList.get(1).equals("+") && !itemList.get(1).equals("-")) {
					throw new RuntimeException("invalid strand! in line:\t" + line);
				}
				if (!IUPACCode.validateNucleotideCode(itemList.get(4))) {
					throw new RuntimeException("invalid character in sequence!\t" + line);
				}
				// h1, i90
				MappedRead mappedRead = new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
						Integer.parseInt(itemList.get(2)),
						itemList.get(4), itemList.get(5));
				mappedRead.generateCpGsInRead(refMap);
				return true;
			}

			@Override
			public Object getResult() {
				return null;
			}
		});

		List<String> resultList = refCpGList.stream()
				.sorted(RefCpG::compareTo)
				.map(r -> String.format("%s\t%d\t%d\t%d\t%f", chr, r.getPos(), r.getCpGCoverage(), r.getMethylCount(),
						r.getMethylLevel()))
				.collect(Collectors.toList());

		Files.asCharSink(new File(outputFileName), Charsets.UTF_8).writeLines(resultList);
	}


}
