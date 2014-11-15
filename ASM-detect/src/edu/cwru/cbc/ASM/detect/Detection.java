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

	public static List<String> execute(String inputName, int threadNumber,
									   Detection detection) throws ExecutionException, InterruptedException {
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
					if (file.isFile() && file.getName().startsWith("chr") && !file.getName().endsWith("aligned") &&
							!file.getName().endsWith("intervalSummary") && !file.getName().endsWith("~") &&
							!file.getName().endsWith(".alignedGroups")) {
						detection.setInputFile(file);
						Future<String> future = executor.submit(detection);
						futureList.add(future);

					}
				}
			}
		} else {
			detection.setInputFile(inputFile);
			Future<String> future = executor.submit(detection);
			futureList.add(future);
		}

		for (Future<String> stringFuture : futureList) {
			resultList.add(stringFuture.get());
		}
		executor.shutdown();

		return resultList;
	}

	public void setInputFile(File inputFile) {
		this.inputFile = inputFile;
	}
}
