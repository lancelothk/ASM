package edu.cwru.cbc.ASM.detect;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * Created by kehu on 11/12/14.
 * Super Class different detection methods.
 */
public abstract class Detection implements Callable<String> {
    protected File inputFile;
    protected String chr;
    protected int startPos;
    protected int endPos;

    public List<String> execute(String inputName, int threadNumber) throws ExecutionException, InterruptedException {
        File inputFile = new File(inputName);
        List<String> resultList = new ArrayList<>();
        ExecutorService executor = Executors.newFixedThreadPool(threadNumber);
        List<Future<String>> futureList = new ArrayList<>();
        if (inputFile.isDirectory()) {
            File[] files = inputFile.listFiles();
            if (files == null) {
                throw new RuntimeException("Empty folder!");
            } else {
                for (File file : files) {
                    try {
                        if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
                                !file.getName().endsWith("test") && !file.getName().endsWith("~") &&
                                !file.getName().endsWith("group")) {
                            Detection newInsDetection = constructNewInstance();
                            newInsDetection.inputFile = file;
                            Future<String> future = executor.submit(newInsDetection);
                            futureList.add(future);
                        }
                    } catch (Exception e) {
                        throw new RuntimeException("Problem File name: " + file.getAbsolutePath() + "\n", e);
                    }
                }
            }
        } else {
            try {
                this.inputFile = inputFile;
                Future<String> future = executor.submit(this);
                futureList.add(future);
            } catch (Exception e) {
                throw new RuntimeException("Problem File name:" + inputFile.getAbsolutePath(), e);
            }
        }

        for (Future<String> stringFuture : futureList) {
            resultList.add(stringFuture.get());
        }
        executor.shutdown();

        return resultList;
    }

    protected void extractIntervalPosition(File inputFile) {
        String[] items = inputFile.getName().split("-");
        if (items.length != 3) {
            throw new RuntimeException("invalid input file name format!\t" + inputFile.getName());
        }
        chr = items[0];
        startPos = Integer.parseInt(items[1]);
        endPos = Integer.parseInt(items[2]);
    }

    abstract protected Detection constructNewInstance();
}
