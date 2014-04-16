package ASM.Execution;

import ASM.DataType.CpGSite;
import ASM.DataType.MappedRead;
import ASM.DataType.MappedReadComparaterByCpG;
import ASM.Utils.MappedReadFileLineProcessor;
import com.google.common.base.Charsets;
import com.google.common.io.Files;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by ke on 2/19/14.
 */
public class FindASM {
    public static void main(String[] args) throws IOException {
        FindASM findASM = new FindASM();
        findASM.execute("/home/ke/ASM_result/chr20/chr20-5727178-5728657.reads", 5727178);
    }

    public void execute(String intervalFileName, long initPos) throws IOException {
        String reference = readRef(intervalFileName);
        List<CpGSite> cpgList = extractCpGSite(reference, initPos);
        List<MappedRead> mappedReadList = Files.asCharSource(new File(intervalFileName), Charsets.UTF_8).readLines(
                new MappedReadFileLineProcessor());
        associateReadWithCpG(cpgList, mappedReadList);
        for (int i = 0; i < mappedReadList.size(); i++) {
            if (mappedReadList.get(i).getCpgList() == null) {
                mappedReadList.remove(i--);
            }
        }
        Collections.sort(mappedReadList, new MappedReadComparaterByCpG());
        findASM(mappedReadList);
    }

    private void findASM(List<MappedRead> mappedReadList) {
        boolean[][] matrix = new boolean[mappedReadList.size()][mappedReadList.size()];
        for (int i = 0; i < mappedReadList.size(); i++) {
            for (int j = 0; j < mappedReadList.size(); j++) {
                if (i != j && comparable(mappedReadList.get(i), mappedReadList.get(j))) {
                    matrix[i][j] = true;
                }
            }
        }
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%d-%d-%s\t", mappedReadList.get(i).getId(), mappedReadList.get(j).getId(),
                                  matrix[i][j]);
            }
            System.out.println();
        }
    }

    private boolean comparable(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            for (CpGSite cpgA : readA.getCpgList()) {
                for (CpGSite cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos() && cpgA.isMethylated() != cpgB.isMethylated()) {
                        return false;
                    }
                }
            }
            return true;
        } else {
            // don't have overlapped CpG, ignore
            return false;
        }
    }

    private String readRef(String intervalFileName) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(intervalFileName));
        String line = bufferedReader.readLine();
        String[] items = line.split("\t");
        if (items.length != 2) {
            throw new RuntimeException("invalid reference line in interval read file!");
        } else {
            return items[1];
        }
    }

    private void associateReadWithCpG(List<CpGSite> cpgList, List<MappedRead> mappedReadList) {
        for (MappedRead read : mappedReadList) {
            for (CpGSite cpg : cpgList) {
                if (cpg.getPos() >= read.getStart() && cpg.getPos() + 1 <= read.getEnd()) {
                    cpg.addAssociatedRead(read);
                    read.addCpG(new CpGSite(cpg.getPos(), read.getCpGMethylStatus(cpg.getPos())));
                }
            }
        }
    }

//    public void generateCpGView(String reference, String intervalFileName, String outputFileName,
//                                       long initPos) throws IOException {
//        List<CpGSite> cpgList = extractCpGSite(reference, initPos);
//        List<MappedRead> mappedReadList = Files.asCharSource(new File(intervalFileName), Charsets.UTF_8).readLines(
//                new MappedReadFileLineProcessor());
//        for (MappedRead mappedRead : mappedReadList) {
//            List<CpGSite> readCpGList = new ArrayList<>();
//            for (CpGSite cpg : cpgList) {
//                if (cpg.getPos() >= mappedRead.getStart() && cpg.getPos() <= mappedRead.getEnd()) {
//                    if (mappedRead.getSequence().charAt((int) (cpg.getPos() - mappedRead.getStart())) == 'C' &&
//                            mappedRead.getSequence().charAt((int) (cpg.getPos() - mappedRead.getStart()) + 1) == 'G') {
//                        readCpGList.add(new CpGSite(cpg.getPos(), true));
//                    }
//                }
//            }
//            mappedRead.setCpgList(readCpGList);
//        }
//        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
//        StringBuilder stringBuilder = new StringBuilder();
//        for (CpGSite cpGSite : cpgList) {
//            stringBuilder.append("\t\t");
//            stringBuilder.append(cpGSite.getPos() - initPos);
//        }
//        bufferedWriter.write(stringBuilder.toString() + "\n");
//        for (MappedRead mappedRead : mappedReadList) {
//            stringBuilder = new StringBuilder();
//            stringBuilder.append(mappedRead.getId());
//            for (CpGSite cpGSite : cpgList) {
//                for (CpGSite mappedCpGSite : mappedRead.getCpgList()) {
//                    if (cpGSite.getPos() == mappedCpGSite.getPos()) {
//                        stringBuilder.append("\t");
//                        stringBuilder.append(mappedCpGSite.isMethylated() ? '*' : '-');
//                    } else {
//                        stringBuilder.append("\t");
//                    }
//                }
//            }
//            bufferedWriter.write(stringBuilder.toString() + "\n");
//        }
//        bufferedWriter.close();
//    }

    private List<CpGSite> extractCpGSite(String reference, long initPos) {
        reference = reference.replace(" ", "");
        List<CpGSite> cpgList = new ArrayList<>();
        for (int i = 0; i < reference.length() - 1; i++) {
            if (reference.charAt(i) == 'c' && reference.charAt(i + 1) == 'g') {
                cpgList.add(new CpGSite(initPos + i));
            }
        }
        return cpgList;
    }
}
