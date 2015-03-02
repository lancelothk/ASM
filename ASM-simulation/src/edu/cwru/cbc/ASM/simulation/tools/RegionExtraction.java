package edu.cwru.cbc.ASM.simulation.tools;

import java.io.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by kehu on 2/16/15.
 * extract reads by regions.
 */
public class RegionExtraction {
    public static void main(String[] args) throws IOException {
        BufferedReader readReader = new BufferedReader(new FileReader("/home/kehu/experiments/ASM/data/i90_r1_chr20"));
        List<String> lines = readReader.lines().collect(Collectors.toList());
        readReader.close();

        BufferedReader regionReader = new BufferedReader(new FileReader(
                "/home/kehu/experiments/ASM/simulation/CpGIslandsRegions/cpgIslandExt_hg18_UCSCGB_chr20_part.txt"));

        String line;
        while ((line = regionReader.readLine()) != null) {
            String[] items = line.split("\t");
            if (items.length <= 3) {
                System.out.println(line);
                continue;
            }
            int start = Integer.parseInt(items[1]);
            int end = Integer.parseInt(items[2]);
            BufferedWriter writer = new BufferedWriter(new FileWriter(
                    String.format("/home/kehu/experiments/ASM/simulation/CpGIslandsRegions/i90_r1_chr20_%d-%d", start,
                                  end)));

            for (String read : lines) {
                items = read.split("\t");
                int x = Integer.parseInt(items[2]);
                int y = Integer.parseInt(items[3]);
                if (x <= end && y >= start) {
                    writer.write(read + "\n");
                }
            }

            writer.close();
        }

        regionReader.close();
    }
}
