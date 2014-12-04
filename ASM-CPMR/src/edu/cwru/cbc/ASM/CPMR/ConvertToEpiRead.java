package edu.cwru.cbc.ASM.CPMR;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import edu.cwru.cbc.ASM.CPMR.DataType.RefChr;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.MappedReadLineProcessorWithSummary;
import edu.cwru.cbc.ASM.commons.DataType.RefCpG;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static edu.cwru.cbc.ASM.commons.Utils.extractCpGSite;

/**
 * Created by lancelothk on 11/14/14.
 * Convert original mapped read into Epiread format
 */
public class ConvertToEpiRead {

    public static void main(String[] args) throws IOException {
        String ref = "hg18_chr20.fa";
        String cellLine = "i90";
        String replicate = "r1";
        String experimentPath = "/home/kehu/experiments/ASM/";

        String referenceGenomeFileName = experimentPath + "/data/ref/" + ref;
        String mappedReadFileName = String.format("%s/data/%s_%s_chr20", experimentPath, cellLine, replicate);
        String outputPath = String.format("%s/result_%s_%s/intervals_chr20", experimentPath, cellLine, replicate);

        convertToEpiRead(referenceGenomeFileName, mappedReadFileName, outputPath);
    }

    public static void convertToEpiRead(String referenceGenomeFileName, String mappedReadFileName,
                                        String outputPath) throws IOException {
        File outputFile = new File(outputPath);
        if (!outputFile.exists()) {
            outputFile.mkdirs();
        }
        RefChr refChr = Utils.readReferenceGenome(referenceGenomeFileName);
        List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), 0);
        List<MappedRead> mappedReadList = Files.readLines(new File(mappedReadFileName), Charsets.UTF_8,
                                                          new MappedReadLineProcessorWithSummary(refCpGList,
                                                                                                 refChr.getRefString().length()));
        Utils.writeExtEpireadInInterval(outputPath, refChr, 0, 0, mappedReadList);
    }
}
