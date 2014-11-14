package edu.cwru.cbc.ASM.detect;

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
        // test reads setting
//        String pathName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/chrTest2-1-6";
//        String summaryFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.summary";
//        String groupResultFileName = "/home/kehu/IdeaProjects/ASM/ASM-detect/testData/test.groupResult";

        // TODO mkdir if not exist
        String cellLine = "i90";
        String replicate = "r1";
        String name = "chr20";

//        String fileName = "chr20-29294521-29294805";
        String fileName = "";

        String pathName = String.format("/home/kehu/experiments/ASM/result_%s_%s/intervals_%s/%s", cellLine, replicate,
                                        name, fileName);

        String summaryFileName = String.format(
                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_%3$s_ASM_summary_%4$s", cellLine, replicate,
                name, fileName);
        String groupResultFileName = String.format(
                "/home/kehu/experiments/ASM/result_%1$s_%2$s/%1$s_%2$s_%3$s_ASM_groups_%4$s", cellLine, replicate, name,
                fileName);

        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(summaryFileName));
        summaryWriter.write("name\tlength\treadCount\tCpGCount\tGroupCount\tavgGroupPerCpG\n");
        BufferedWriter groupWriter = new BufferedWriter(new FileWriter(groupResultFileName));
        File path = new File(pathName);

        ExecutorService executor = Executors.newFixedThreadPool(1);
        List<Future<String>> futureList = new ArrayList<>();
        if (path.isDirectory()) {
            for (File file : path.listFiles()) {
                if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
                        !file.getName().endsWith("intervalSummary") && !file.getName().endsWith("~") &&
                        !file.getName().endsWith(".alignedGroups")) {
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
//            System.out.println(stringFuture.get());
        }
        executor.shutdown();

        summaryWriter.close();
        groupWriter.close();
        System.out.println(System.currentTimeMillis() - start + "ms");
    }
}