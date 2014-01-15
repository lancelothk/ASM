import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class Alignment {
	public static void main(String[] args) {
		
	}
	
	public static ArrayList<Read> readMappedReads(String fileName){
		ArrayList<Read> readsList = new ArrayList<>();
		try {
			BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
			String line=null;
			while ((line =bufferedReader.readLine()) != null){
				String[] items = line.split("\t");
				readsList.add(new Read(items[0],items[1],Long.parseLong(items[2]), Long.parseLong(items[3]),items[4]));
			}
			bufferedReader.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return readsList;
	}
	
	public static void alignReads(List<Read> readsList){
		// sort reads first
		Collections.sort(readsList, new ReadComparator());
		
	}
}
