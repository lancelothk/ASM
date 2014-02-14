package ASM;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Alignment {
	public static void main(String[] args) {
		String targetFileName = "/media/ke/win-data/Dataset/WholeGenomeMethylation/CPGI_chr6_cpg35";
		String outputFileName = targetFileName + ".aligned";
		List<Read> readsList = readMappedReads(targetFileName);
		alignReads(readsList, outputFileName);
	}

	public static ArrayList<Read> readMappedReads(String inputFileName) {
		ArrayList<Read> readsList = new ArrayList<>();
		BufferedReader bufferedReader;
		try {
			bufferedReader = new BufferedReader(new FileReader(inputFileName));
            String line;
            while ((line = bufferedReader.readLine()) != null) {
				String[] items = line.split("\t");
				if (!items[1].equals("+") && !items[1].equals("-")){
					System.err.println("invalid strand symbol!");
				}
				readsList.add(new Read(items[0], items[1].charAt(0), Long.parseLong(items[2]), Long.parseLong(items[3]), items[4], items[5]));
			}
			bufferedReader.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return readsList;
	}

	public static void alignReads(List<Read> readsList, String outputFileName) {
		// sort reads first
		Collections.sort(readsList, new ReadComparator());
		// set initial position
		long initialPos = readsList.get(0).getStart();
		BufferedWriter bufferedWriter;
		try {
			bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
			for (Read read : readsList) {
				bufferedWriter.write(read.toString(initialPos) + "\n");
			}
			bufferedWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
}
