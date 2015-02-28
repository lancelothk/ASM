package edu.cwru.cbc.ASM.tools.Conversion;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.DataType.MappedRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kehu on 2/6/15.
 * convert amrfinder MR format to Lister 2009 mapped reads format.
 */
public class MRToMappedRead {
    public static void main(String[] args) throws IOException {
        String input = "/home/kehu/experiments/ASM/amrfinder/simulation/i90_r1_chr20_sim.mr";
        String output = "/home/kehu/experiments/ASM/amrfinder/simulation/i90_r1_chr20_sim";

        List<MappedRead> mappedReadList = Files.readLines(new File(input), Charsets.UTF_8,
                                                          new LineProcessor<List<MappedRead>>() {
                                                              private List<MappedRead> mappedReadList = new ArrayList<>();

                                                              @Override
                                                              public boolean processLine(
                                                                      String line) throws IOException {
                                                                  String[] items = line.split("\t");
                                                                  if (items.length != 8) {
                                                                      throw new RuntimeException(
                                                                              "invalid mapped read format:" + line);
                                                                  }
                                                                  if (items[5].length() != 1) {
                                                                      throw new RuntimeException("invalid strand!");
                                                                  }

                                                                  MappedRead mappedRead = new MappedRead(items[0],
                                                                                                         items[5].charAt(
                                                                                                                 0),
                                                                                                         Integer.parseInt(
                                                                                                                 items[1]),
                                                                                                         Integer.parseInt(
                                                                                                                 items[2]),
                                                                                                         items[6],
                                                                                                         items[3]);
                                                                  mappedReadList.add(mappedRead);
                                                                  return true;
                                                              }

                                                              @Override
                                                              public List<MappedRead> getResult() {
                                                                  return mappedReadList;
                                                              }
                                                          });

        BufferedWriter mappedReadWriter = new BufferedWriter(new FileWriter(output));
        for (MappedRead mappedRead : mappedReadList) {
            mappedReadWriter.write(mappedRead.outputString());
        }
        mappedReadWriter.close();
    }
}
