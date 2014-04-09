package ASM.DataType;

import ASM.Utils.Utils;
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
    private String chr;
    private String intervalFolderName;

    public IntervalMatchingLineProcessor(ChrCoverageSummary chrCoverageSummary, String intervalFolderName,
                                         String chr) throws IOException {
        System.out.println("start interval matching:");
        this.intervalList = chrCoverageSummary.generateIntervals();
        this.chr = chr;
        this.intervalFolderName = intervalFolderName;
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
                    interval.addRead(line);
                }
            }
            counter++;
            if (counter % 1000000 == 0) {
                System.out.printf("%d\t%s%n", counter, Utils.getCurrentTime());
            }
            return true;
        }
    }

    @Override
    public List<GenomicInterval> getResult() {
        System.out.println(String.format("Average Read Length:%f\nMax Read Length:%d\nMin Read Length:%d",
                                         totalReadLength / (double) totalReadCount, maxReadLength, minReadLength));
        File intervalFolder = new File(intervalFolderName);
        if (!intervalFolder.exists()) {
            intervalFolder.mkdir();
        }
        System.out.println("start to write reads for intervals\t" + Utils.getCurrentTime());
        for (GenomicInterval interval : intervalList) {
            try {
                Utils.writeReads(interval.getReadList(),
                                 String.format("%s/%s-%d-%d.reads", intervalFolderName, chr, interval.getStart(),
                                               interval.getEnd())
                                );
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return this.intervalList;
    }
}
