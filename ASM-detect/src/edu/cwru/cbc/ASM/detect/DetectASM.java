package edu.cwru.cbc.ASM.detect;

import edu.cwru.cbc.ASM.align.AlignReads;
import edu.cwru.cbc.ASM.detect.WithMappedRead.DetectionWihMappedRead;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Created by ke on 2/19/14.
 * ASM detection
 * Note:    1. mapped reads start pos is 1-based, end pos is 0-based.
 */
public class DetectASM {

    public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
        long start = System.currentTimeMillis();
        // call on single file
//        detectASM.execute("/home/kehu/ASM_result/chr20-56897421-56898208.reads", 56897421);
//        detectASM.execute("", 1);
//		detectASM.execute(new File("/home/lancelothk/chr20_test/chr20-56895353-56895567"), 56895353);
        // test reads setting
//        String pathName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/chrTest2-1-6";
//        String summaryFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.summary";
//        String groupResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.groupResult";
//        String group2ResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.group2Result";

        // TODO mkdir if not exist
        String cellLine = "i90";
        String replicate = "r1";
        String name = "_chr20";

        String pathName = String.format("/home/kehu/experiments/ASM/result_%s_%s/intervals%s/chr20-56850168-56850901",
                                        cellLine, replicate, name);
        AlignReads.align(pathName);

        String summaryFileName = String.format(
                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_summary%3$s", cellLine, replicate,
                name);
        String groupResultFileName = String.format(
                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_groups%3$s", cellLine, replicate,
                name);
        String group2ResultFileName = String.format(
                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_chr22_ASM_group2%3$s", cellLine, replicate,
                name);

        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
        summaryWriter.write("name\tlength\treadCount\tCpGCount\tGroupCount\tavgGroupPerCpG\n");
        BufferedWriter groupWriter = new BufferedWriter(new FileWriter(groupResultFileName));
        File path = new File(pathName);

        ExecutorService executor = Executors.newFixedThreadPool(4);
        List<Future<String>> futureList = new ArrayList<>();
        if (path.isDirectory()) {
            for (File file : path.listFiles()) {
                if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
                        !file.getName().endsWith("intervalSummary") && !file.getName().endsWith("~")) {
                    String[] items = file.getName().split("-");
                    Future<String> future = executor.submit(
                            new DetectionWihMappedRead(file, Integer.parseInt(items[1]), groupWriter));
                    futureList.add(future);
                }
            }
        } else {
            String[] items = path.getName().split("-");
            Future<String> future = executor.submit(
                    new DetectionWihMappedRead(path, Integer.parseInt(items[1]), groupWriter));
            futureList.add(future);
        }

        for (Future<String> stringFuture : futureList) {
            summaryWriter.write(stringFuture.get());
            //System.out.println(stringFuture.get());
        }
        executor.shutdown();

        summaryWriter.close();
        groupWriter.close();
        System.out.println(System.currentTimeMillis() - start + "ms");
    }
}