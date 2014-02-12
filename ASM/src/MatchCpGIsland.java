import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MatchCpGIsland {

	public static void main(String[] args) throws IOException {
		MatchCpGIsland.matchMappedReadToCpGIsland("/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/reads_bs_i90_r1.mapped_chr6",
				"/media/ke/win-data/Dataset/cpgIslandExt_hg18_UCSCGB_chr6.txt", "/home/ke/test");
	}
	
	public static void matchMappedReadToCpGIsland(String mappedReadFileName, String cpgIslandFileName, String outputPath) throws IOException{
		List<CpGIsland> cpGIslandsList = readCpGIslands(cpgIslandFileName);
		List<MappedRead> mappedReadsList = readMappedReads(mappedReadFileName);
		for (CpGIsland cpGIsland : cpGIslandsList) {
			BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputPath + "/" + cpGIsland.getName()));
			for (MappedRead mappedRead : mappedReadsList) {
				if (readInIsland(mappedRead, cpGIsland)){
					bufferedWriter.write(mappedRead.toRead());
				}
			}
			bufferedWriter.close();
			System.out.println(cpGIsland.getName() + "\t finished");
		}
	}

	public static List<CpGIsland> readCpGIslands(String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		List<CpGIsland> cpgIslandList = new ArrayList<>();
		String line = null;
		String[] items;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split("\t");
			cpgIslandList.add(new CpGIsland(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]), items[0] + items[1], Integer
					.parseInt(items[4]), Integer.parseInt(items[5]), Integer.parseInt(items[6]), Double.parseDouble(items[7]), Double
					.parseDouble(items[8]), Double.parseDouble(items[9])));
		}
		bufferedReader.close();
		System.out.println("readCpGIslands finished");
		return cpgIslandList;
	}

	public static List<MappedRead> readMappedReads(String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
		List<MappedRead> mappedReadsList = new ArrayList<>();
		String line = null;
		String[] items;
		long counter = 0;
		while ((line = bufferedReader.readLine()) != null) {
			items = line.split("\t");
			if (items[0].startsWith("chrom")){
				continue;
			}
			mappedReadsList.add(new MappedRead(items[0], items[1], Long.parseLong(items[2]), Long.parseLong(items[3]), items[4], Long
					.parseLong(items[5])));
			counter++;
			if (counter %100000 == 0){
				System.out.println(counter);
			}
		}
		bufferedReader.close();
		System.out.println("readMappedReads finished");
		return mappedReadsList;
	}
	
	public static boolean readInIsland(MappedRead read, CpGIsland island){
		if ((read.getStart() <= island.getEnd() && read.getStart() >= island.getStart()) || (read.getEnd() <= island.getEnd() && read.getEnd() >= island.getStart())){
			return true;
		}else {
			return false;
		}
	}
}
