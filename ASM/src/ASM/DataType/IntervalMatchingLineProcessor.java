package ASM.DataType;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.io.LineProcessor;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by ke on 3/5/14.
 * Implement interval and read matching
 */
public class IntervalMatchingLineProcessor implements LineProcessor<List<GenomicInterval>> {
    private List<GenomicInterval> intervalList;
    private int counter = 0;
    private long totalReadLength;
    private long totalReadCount;
    private int minReadLength = Integer.MAX_VALUE;
    private int maxReadLength = Integer.MIN_VALUE;

    public IntervalMatchingLineProcessor(ChrCoverageSummary chrCoverageSummary, String intervalFolerName,
                                         String chr) throws IOException {
        System.out.println("start interval matching:");
        this.intervalList = chrCoverageSummary.generateIntervals();
        File intervalFolder = new File(intervalFolerName);
        if (!intervalFolder.exists()) {
            intervalFolder.mkdir();
        }
        for (GenomicInterval genomicInterval : intervalList) {
            genomicInterval.initializeWriter(intervalFolerName, chr);
        }
    }

    @Override
    public boolean processLine(String line) throws IOException {
        if (line.startsWith("chr")) {
            return true;
        } else if (line.equals("")) {
            return false;
        } else {
            List<String> itemList = Lists.newArrayList(Splitter.on('\t').split(line));
            int start = Integer.parseInt(itemList.get(2));
            int end = Integer.parseInt(itemList.get(3));
            int length = end - start + 1;
            totalReadLength += length;
            totalReadCount++;
            if (length > maxReadLength) {
                maxReadLength = length;
            }
            if (length < minReadLength) {
                minReadLength = length;
            }

            for (GenomicInterval interval : intervalList) {
                if (start >= interval.getStart() && interval.getEnd() >= end) {
                    interval.incrementReadCount();
                    interval.write(line + "\n");
                }
            }
            counter++;
            if (counter % 1000000 == 0) {
                System.out.println(counter);
            }
            return true;
        }
    }

    @Override
    public List<GenomicInterval> getResult() {
        System.out.println(String.format("Average Read Length:%f\nMax Read Length:%d\nMin Read Length:%d",
                                         totalReadLength / (double) totalReadCount, maxReadLength, minReadLength));
        for (GenomicInterval genomicInterval : intervalList) {
            try {
                genomicInterval.closeWriter();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return this.intervalList;
    }
}
