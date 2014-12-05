package edu.cwru.cbc.ASM.commons;

import java.io.*;
import java.util.logging.Logger;

/**
 * Created by ke on 2/19/14.
 * Split original downloaded whole file into chromosomes.
 */
public class SplitFileByChr {
    private static final Logger logger = Logger.getLogger(SplitFileByChr.class.getName());

    public static void main(String[] args) throws IOException {
        String inputFileName = "/media/ke/win-data/Dataset/WholeGenomeMethylation/reads_bs_i90_r1.mapped";
        String outputFilePath = "/media/ke/win-data/Dataset/WholeGenomeMethylation/i90_r1/";
        split(inputFileName, outputFilePath);
    }

    public static void split(String inputFileName, String outputFilePath) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFileName));
        BufferedWriter bufferedWriter = null;
        String line, chr = "";
        String[] items;
        while ((line = bufferedReader.readLine()) != null) {
            if (line.startsWith("chrom")) {
                continue;
            }
            items = line.split("\t");
            if (!chr.equals(items[0])) {
                if (bufferedWriter != null) {
                    logger.info(chr + "\tfinished");
                    bufferedWriter.close();
                }
                chr = items[0];
                bufferedWriter = new BufferedWriter(new FileWriter(String.format("%schr%s", outputFilePath, chr)));
                bufferedWriter.write(String.format("%s\n", line));
            } else {
                assert bufferedWriter != null;
                bufferedWriter.write(String.format("%s\n", line));
            }
        }
        if (bufferedWriter != null) {
            bufferedWriter.close();
        }
        bufferedReader.close();
    }
}
