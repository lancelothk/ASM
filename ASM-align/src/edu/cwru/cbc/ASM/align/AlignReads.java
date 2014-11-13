package edu.cwru.cbc.ASM.align;

import edu.cwru.cbc.ASM.commons.DataType.MappedRead;
import edu.cwru.cbc.ASM.commons.DataType.ReadComparator;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Align reads to readable format with padding.
 */

public class AlignReads {
    private static String ref;

    public static void main(String[] args) throws Exception {
        String targetFileName;
        if (args.length == 0) {
            targetFileName = "/home/kehu/experiments/ASM/result_i90_r1/intervals_newParameter_large/chr22-15240456-15242527";
        } else {
            targetFileName = args[0];
        }
        File targetFile = new File(targetFileName);
        if (targetFile.isDirectory()) {
            File[] files = targetFile.listFiles();
            assert files != null;
            for (File file : files) {
                try {
                    if (file.isFile() && !file.isHidden() && !file.getName().endsWith(".aligned")) {
                        String outputFileName = file.getAbsolutePath() + ".aligned";
                        List<MappedRead> readsList = readMappedReads(file);
                        alignReads(readsList, outputFileName);
                    }
                } catch (Exception e) {
                    throw new Exception(file.getAbsolutePath(), e);
                }
            }
        } else {
            String outputFileName = targetFileName + ".aligned";
            List<MappedRead> readsList = readMappedReads(targetFile);
            alignReads(readsList, outputFileName);
        }
    }

    public static ArrayList<MappedRead> readMappedReads(File inputFile) {
        ArrayList<MappedRead> readsList = new ArrayList<>();
        BufferedReader bufferedReader;
        try {
            bufferedReader = new BufferedReader(new FileReader(inputFile));
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                if (line.startsWith("ref")) {
                    ref = line.split("\t")[1];
                    continue;
                }
                String[] items = line.split("\t");
                if (!items[1].equals("+") && !items[1].equals("-")) {
                    System.err.println("invalid strand symbol!");
                }
                readsList.add(new MappedRead(items[0], items[1].charAt(0), Integer.parseInt(items[2]),
                                             Integer.parseInt(items[3]),
                                       items[4], items[5]));
            }
            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }
        return readsList;
    }

    public static void alignReads(List<MappedRead> readsList, String outputFileName) {
        // sort reads first
        readsList.sort(new ReadComparator());
        // set initial position
        int initialPos = readsList.get(0).getStart();
        BufferedWriter bufferedWriter;
        try {
            bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
            if (ref != null) {
                bufferedWriter.write(String.format("ref:\t%s\n", ref));
            }
            for (MappedRead read : readsList) {
                bufferedWriter.write(read.toString(initialPos) + "\n");
            }
            bufferedWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
