package ASM.Execution;

import ASM.DataType.CpGSite;
import ASM.DataType.MappedRead;
import ASM.Utils.IOUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/19/14.
 */
public class FindASM {
	public static void main(String[] args) throws IOException {
		final String intervalRef = "gcgtggccgcggggcccgtgctggctgggaccacctccaggagcccgtgg gtgcccgggaagagggcaggaggtgctgggccctgctggggaagaggatc ggtaaggacttgcaggtgactttgctcgcctgtcctgaacagggcacagg gaaggcttcccagggccgtctccttgtggcctggctcctgcgggccccga gccagcccttcccggccagcccagcacaacagccccgccagtccgaacca atggagaactgcttcgcgatgggaagccgacaccagcagcgaaagcacaa cttccggccaggggaaaagcacacctccaaccatgcggcgctaggcctcg gctctccagtccgcaaggtcctttccaggccgccgcgggcgctgctgagc tctggaccggactgtcccgccctcccc";
		final String readsFileName = "/media/ke/win-data/Dataset/WholeGenomeMethylation/ASM/CPGI_chr6_cpg35";
		final String outputFileName = "/home/ke/test/out";
		final long initPos = 368883;
		generateCpGView(intervalRef,readsFileName,outputFileName, initPos);
	}

	public static void find(){
		// TODO implement find method
	}

	public static void generateCpGView(String reference, String intervalFileName, String outputFileName, long initPos) throws IOException {
		List<CpGSite> cpgList = extractCpGSite(reference, initPos);
		List<MappedRead> mappedReadList = IOUtils.readMappedRead(intervalFileName);
		for (MappedRead mappedRead : mappedReadList) {
			List<CpGSite> readCpGList = new ArrayList<>();
			for (CpGSite cpg : cpgList) {
				if (cpg.getPos() >= mappedRead.getStart() && cpg.getPos() <= mappedRead.getEnd()){
					if (mappedRead.getSequence().charAt((int) (cpg.getPos() - mappedRead.getStart())) == 'C'
							&& mappedRead.getSequence().charAt((int) (cpg.getPos() - mappedRead.getStart()) + 1) == 'G'){
						readCpGList.add(new CpGSite(cpg.getPos(), true));
					}
				}
			}
			mappedRead.setCpgList(readCpGList);
		}
		BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputFileName));
		StringBuilder stringBuilder = new StringBuilder();
		for (CpGSite cpGSite : cpgList) {
			stringBuilder.append("\t\t");
			stringBuilder.append(cpGSite.getPos() - initPos);
		}
		bufferedWriter.write(stringBuilder.toString() + "\n");
		for (MappedRead mappedRead : mappedReadList) {
			stringBuilder = new StringBuilder();
			stringBuilder.append(mappedRead.getId());
			for (CpGSite cpGSite : cpgList) {
				for (CpGSite mappedCpGSite : mappedRead.getCpgList()) {
					if (cpGSite.getPos() == mappedCpGSite.getPos()){
						stringBuilder.append("\t");
						stringBuilder.append(mappedCpGSite.isMethylated()?'*':'-');
					}else {
						stringBuilder.append("\t");
					}
				}
			}
			bufferedWriter.write(stringBuilder.toString() + "\n");
		}
		bufferedWriter.close();
	}

	private static List<CpGSite> extractCpGSite(String reference, long initPos){
		reference = reference.replace(" ", "");
		List<CpGSite> cpgList = new ArrayList<>();
		for (int i = 0; i < reference.length()-1; i++) {
			if (reference.charAt(i) == 'c' && reference.charAt(i+1) == 'g'){
				cpgList.add(new CpGSite(initPos + i));
			}
		}
		return cpgList;
	}
}
