package ASM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MatchCpGIsland {

	public static void main(String[] args) throws IOException {
		MatchCpGIsland.matchMappedReadToCpGIsland("/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/reads_bs_i90_r1.mapped_chr6",
				"/media/ke/win-data/Dataset/cpgIslandExt_hg18_UCSCGB_chr6.txt", "/home/ke/test");
	}

	public static void matchMappedReadToCpGIsland(String mappedReadFileName, String cpgIslandFileName, String outputPath) throws IOException {
		List<CpGIsland> cpGIslandsList = readCpGIslands(cpgIslandFileName);
		for (CpGIsland cpGIsland : cpGIslandsList) {
			cpGIsland.initWriter(outputPath);
		}
		BufferedReader bufferedReader = new BufferedReader(new FileReader(mappedReadFileName));
        String line;
        String[] items;
		long counter = 0;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split("\t");
			if (items[0].startsWith("chrom")) {
				continue;
			}
			MappedRead mappedRead = new MappedRead(items[0], items[1], Long.parseLong(items[2]), Long.parseLong(items[3]), items[4],
					Long.parseLong(items[5]));
			for (CpGIsland cpGIsland : cpGIslandsList) {
				if (readInIsland(mappedRead, cpGIsland)) {
					cpGIsland.write(mappedRead.toRead());
				}
			}
			counter++;
			if (counter % 1000000 == 0) {
				System.out.println(counter);
			}
		}
		bufferedReader.close();
		for (CpGIsland cpGIsland : cpGIslandsList) {
			cpGIsland.closeWriter();
		}
	}

	public static List<CpGIsland> readCpGIslands(String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		List<CpGIsland> cpgIslandList = new ArrayList<>();
        String line;
        String[] items;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split("\t");
			cpgIslandList.add(new CpGIsland(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]), String.format("%s-%s", items[0], items[1]), Integer
					.parseInt(items[4]), Integer.parseInt(items[5]), Integer.parseInt(items[6]), Double.parseDouble(items[7]), Double
					.parseDouble(items[8]), Double.parseDouble(items[9])));
		}
		bufferedReader.close();
		System.out.println("readCpGIslands finished");
		return cpgIslandList;
	}

	public static boolean readInIsland(MappedRead read, CpGIsland island) {
        return (read.getStart() <= island.getEnd() && read.getStart() >= island.getStart())
                || (read.getEnd() <= island.getEnd() && read.getEnd() >= island.getStart());
    }
}
