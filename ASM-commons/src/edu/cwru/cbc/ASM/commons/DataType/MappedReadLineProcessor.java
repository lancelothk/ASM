package edu.cwru.cbc.ASM.commons.DataType;

import com.google.common.io.LineProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by ke on 2/21/14.
 * Implemented LineProcessor for reading MappedRead File
 * Mapped reads start pos is 1-based, end pos is 0-based.
 */
public class MappedReadLineProcessor implements LineProcessor<List<MappedRead>> {
    private static final int MIN_READ_CPG = 2;
    protected List<MappedRead> mappedReadList = new ArrayList<>();
    protected Map<Integer, RefCpG> refMap;

    public MappedReadLineProcessor(List<RefCpG> refCpGList) {
        this.refMap = refCpGList.stream().collect(Collectors.toMap(RefCpG::getPos, refCpG -> refCpG));
    }

    @Override
    public boolean processLine(String line) throws IOException {
        try {
            if (line.startsWith("chr") || line.startsWith("ref")) {
                return true;
            } else if (line.equals("")) {
                return false;
            } else {
                MappedRead mappedRead = processRead(line);
                if (mappedRead.getCpgList().size() >= MIN_READ_CPG) {
                    mappedReadList.add(mappedRead);
                }
                return true;
            }
        } catch (Exception e) {
            throw new RuntimeException("Problem line: " + line + "\n", e);
        }
    }

    protected MappedRead processRead(String line) {
        String[] items = line.split("\t");
        if (items.length != 6) {
            throw new RuntimeException("invalid mapped read format:" + line);
        }
        if (items[1].length() != 1) {
            throw new RuntimeException("invalid strand!");
        }
        int start = Integer.parseInt(items[2]);// mapped read is 1bp right than UCSC ref
        int end = Integer.parseInt(items[3]);

        MappedRead mappedRead = new MappedRead(items[0], items[1].charAt(0), Integer.parseInt(items[2]),
                                               Integer.parseInt(items[3]), items[4], items[5]);
        for (int i = start; i < end; i++) {
            if (refMap.containsKey(i)) {
                CpG cpg = new CpG(mappedRead, refMap.get(i), mappedRead.getMethylStatus(i));
                refMap.get(i).addCpG(cpg);
                // ignore unknown methyl CpG
//  				if (cpg.getMethylStatus() != MethylStatusq.N) {
                mappedRead.addCpG(cpg);
//	    			}
                i++;
            }
        }
        return mappedRead;
    }

    @Override
    public List<MappedRead> getResult() {
        return mappedReadList;
    }
}
