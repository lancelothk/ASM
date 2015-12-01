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
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 * Created by kehu on 12/1/15.
 * Merge CpG from +/- strand into one. The input is from bismark coverage output.
 */
public class MergeDiffStrandCpGPgm {
	public static void main(String[] args) throws IOException, ParseException, ExecutionException,
			InterruptedException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input cov file").required().build());
		options.addOption(Option.builder("o").hasArg().desc("output file").required().build());
		options.addOption(Option.builder("r").hasArg().desc("Reference File path").required().build());

		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);

		String inputFileName = cmd.getOptionValue("i");
		String referenceFileName = cmd.getOptionValue("r");
		String outputFileName = cmd.getOptionValue("o");
		execution(inputFileName, referenceFileName, outputFileName);
	}

	private static void execution(String inputFileName, String referenceFileName, String outputFileName) throws
			IOException {
		Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap = MethylationUtils.initializeGenomeRefCpGMap(
				IOUtils.readReferenceGenome(referenceFileName));
		Files.readLines(new File(inputFileName), Charsets.UTF_8, new LineProcessor<Object>() {
			@Override
			public boolean processLine(@Nonnull String line) throws IOException {
				return false;
			}

			@Override
			public Object getResult() {
				return null;
			}
		});


	}
}
