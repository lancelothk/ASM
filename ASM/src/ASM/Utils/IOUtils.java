package ASM.Utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Created by ke on 4/9/14.
 */
public class IOUtils {

    public static void writeReads(List<String> readList, String fileName) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
        for (String s : readList) {
            writer.write(s + "\n");
        }
        writer.close();
    }
}
