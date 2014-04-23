package ASM.Execution;

import ASM.DataType.CpGSite;
import ASM.DataType.MappedRead;
import ASM.DataType.MappedReadComparaterByCpG;
import ASM.DataType.Node;
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
    private enum Compatibility {
        COMPATIBLE, NONCOMPATIBLE, NONCONNECTED, SELF
    }

    private Compatibility[][] matrix;

    public static void main(String[] args) throws IOException {
        FindASM findASM = new FindASM();
        findASM.execute("/home/ke/ASM_result/chr20/chr20-56861267-56862291.reads", 56861267);
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
        constructGraph(mappedReadList);
        detectASM(mappedReadList);
    }

    private List<List<Node>> detectASM(List<MappedRead> mappedReadList) {
        List<List<Node>> resultList = new ArrayList<>();
        List<Node> nodeList = new ArrayList<>();
        for (int i = 0; i < mappedReadList.size(); i++) {
            nodeList.add(new Node(i));
        }
        while (true) {
            List<Node> result = getGroup(nodeList);
            if (result.size() == 0) {
                break;
            } else {
                resultList.add(result);
                nodeList.removeAll(result);
            }
        }
        return resultList;
    }

    private List<Node> getGroup(List<Node> nodeList) {
        List<Node> result = new ArrayList<>();
        List<Node> seedList = new ArrayList<>();
        seedList.add(nodeList.get(0));
        result.add(nodeList.get(0));
        while (true) {
            if (seedList.size() == 0) {
                break;
            } else {
                addCandidates(nodeList, result, seedList, seedList.get(0));
            }
        }
        return result;
    }

    private void addCandidates(List<Node> nodeList, List<Node> result, List<Node> seedList, Node seed) {
        // check all other reads
        for (int i = 0; i < matrix[seed.getIndex()].length; i++) {
            if (matrix[seed.getIndex()][i] == Compatibility.COMPATIBLE) {
                //if candidate has no contradiction with result member, add to result
                if (!hasContradiction(result, i)) {
                    seedList.add(nodeList.get(i));
                    result.add(nodeList.get(i));
                }
            }
        }
        seedList.remove(seed);
    }

    private boolean hasContradiction(List<Node> result, int candidate) {
        for (Node member : result) {
            if (matrix[candidate][member.getIndex()] == Compatibility.NONCOMPATIBLE) {
                return true;
            }
        }
        return false;
    }

    private void constructGraph(List<MappedRead> mappedReadList) {
        matrix = new Compatibility[mappedReadList.size()][mappedReadList.size()];
        for (int i = 0; i < mappedReadList.size(); i++) {
            for (int j = 0; j < mappedReadList.size(); j++) {
                if (i == j) {
                    matrix[i][j] = Compatibility.SELF;
                }else {
                    matrix[i][j] = checkCompatible(mappedReadList.get(i), mappedReadList.get(j));
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

    private Compatibility checkCompatible(MappedRead readA, MappedRead readB) {
        if (readA.getFirstCpG().getPos() <= readB.getLastCpG().getPos() &&
                readA.getLastCpG().getPos() >= readB.getFirstCpG().getPos()) {
            for (CpGSite cpgA : readA.getCpgList()) {
                for (CpGSite cpgB : readB.getCpgList()) {
                    if (cpgA.getPos() == cpgB.getPos() && cpgA.isMethylated() != cpgB.isMethylated()) {
                        return Compatibility.NONCOMPATIBLE;
                    }
                }
            }
            return Compatibility.COMPATIBLE;
        } else {
            // don't have overlapped CpG, ignore
            return Compatibility.NONCONNECTED;
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
