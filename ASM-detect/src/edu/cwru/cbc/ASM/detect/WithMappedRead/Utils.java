package edu.cwru.cbc.ASM.detect.WithMappedRead;

import org.apache.commons.lang3.StringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by kehu on 12/9/14.
 */
public class Utils {
    public static String readRefFromIntervalFile(File inputFile) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(inputFile));
        String line = bufferedReader.readLine();
        String[] items = line.split("\t");
        if (items.length != 2) {
            throw new RuntimeException(
                    "invalid reference line in interval read file!\t" + line + "\t" + inputFile.getName());
        } else {
            // item[1](reference string) should be in uppercase and without space
            assert StringUtils.isAllUpperCase(items[1]);
            assert !items[1].contains(" ");
            return items[1];
        }
    }
}
