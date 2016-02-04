package edu.cwru.cbc.ASM.tools;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.methylation.CpG;
import edu.cwru.cbc.ASM.commons.methylation.MethylStatus;
import edu.cwru.cbc.ASM.commons.methylation.MethylationUtils;
import edu.cwru.cbc.ASM.commons.methylation.RefCpG;
import edu.cwru.cbc.ASM.commons.sequence.MappedRead;
import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;
import org.apache.commons.cli.*;

import javax.annotation.Nonnull;
import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by lancelothk on 12/23/15.
 * Visualize methylation pattern and SNP
 */
public class MethylFigurePgm {
	private static final Splitter tabSplitter = Splitter.on("\t");
	private static final int COMMON_FONT_SIZE = 16;
	private static final int CG_RADIUS = 20;
	private static final int HEIGHT_INTERVAL = 24;
	private static final int BPWIDTH = 10;

	public static void main(String[] args) throws ParseException, IOException {
		Options options = new Options();
		options.addOption(Option.builder("i").hasArg().desc("input grouped read file").build());
		options.addOption(Option.builder("p").hasArg().desc("SNP position").build());
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = parser.parse(options, args);
		String groupedReadFile = cmd.getOptionValue("i");
		int snpPosition = Integer.parseInt(cmd.getOptionValue("p"));

		Files.readLines(new File(groupedReadFile), Charsets.UTF_8,
				new LineProcessor() {
					String ref;
					int group = 0;
					List<MappedRead> group1 = new ArrayList<>();
					List<MappedRead> group2 = new ArrayList<>();

					@Override
					public boolean processLine(@Nonnull String line) throws IOException {
						List<String> itemList = tabSplitter.splitToList(line);
						if (line.startsWith("ref:")) {
							ref = itemList.get(1);
						} else if (line.equals("")) {
							group++;
						} else {
							// TODO add input data validation
							switch (group) {
								case 0:
									group1.add(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
											Integer.parseInt(itemList.get(2)),
											itemList.get(1).charAt(0) == '+' ? itemList.get(4)
													.replace(".", "") : MappedRead.getComplementarySequence(
													itemList.get(4).replace(".", "")), itemList.get(5)));
									break;
								case 1:
									group2.add(new MappedRead(itemList.get(0), itemList.get(1).charAt(0),
											Integer.parseInt(itemList.get(2)),
											itemList.get(1).charAt(0) == '+' ? itemList.get(4)
													.replace(".", "") : MappedRead.getComplementarySequence(
													itemList.get(4).replace(".", "")), itemList.get(5)));
									break;
								default:
									throw new RuntimeException("more than 2 groups!");
							}
						}
						return true;
					}

					@Override
					public Object getResult() {
						if (group2.size() == 0) {
							int minStart = group1.stream().mapToInt(MappedRead::getStart).min().getAsInt();
							List<RefCpG> refCpGList = MethylationUtils.extractCpGSite(ref, minStart);
							HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
							for (RefCpG refCpG : refCpGList) {
								refMap.put(refCpG.getPos(), refCpG);
							}
							for (MappedRead mappedRead : group1) {
								mappedRead.generateCpGsInRead(refMap);
							}
							drawCompactFigure(group1, snpPosition, refCpGList, groupedReadFile + ".compact");
							return null;
						} else {
							int minStart = Math.min(group1.stream().mapToInt(MappedRead::getStart).min().getAsInt(),
									group2.stream().mapToInt(MappedRead::getStart).min().getAsInt());
							int maxEnd = Math.max(group1.stream().mapToInt(MappedRead::getEnd).max().getAsInt(),
									group2.stream().mapToInt(MappedRead::getEnd).max().getAsInt());
							List<RefCpG> refCpGList = MethylationUtils.extractCpGSite(ref, minStart);
							HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
							for (RefCpG refCpG : refCpGList) {
								refMap.put(refCpG.getPos(), refCpG);
							}
							for (MappedRead mappedRead : group1) {
								mappedRead.generateCpGsInRead(refMap);
							}
							for (MappedRead mappedRead : group2) {
								mappedRead.generateCpGsInRead(refMap);
							}
							drawCompactFigure(group1, group2, snpPosition, groupedReadFile + ".compact");
//						drawFigure(group1, group2, minStart, maxEnd, snpPosition, groupedReadFile + ".png");
							return null;
						}
					}
				});
	}

	private static void drawCompactFigure(List<MappedRead> group1, int snpPosition, List<RefCpG> refCpGList,
	                                      String outputFileName) {
		int imageHeight = (group1.size() + 1) * HEIGHT_INTERVAL, imageWidth = (refCpGList.size() + 3) * CG_RADIUS;

		try {
			EPSWriter epsWriter = new EPSWriter(outputFileName, imageWidth, imageHeight);
			Graphics2D graphWriter = epsWriter.getGraphWriter();
			graphWriter.setBackground(Color.WHITE);
			graphWriter.clearRect(0, 0, imageWidth, imageHeight);
			graphWriter.setPaint(Color.BLACK);
			graphWriter.setFont(new Font("Helvetica", Font.PLAIN, COMMON_FONT_SIZE));

			int height = 0;
			drawCompactGroup(graphWriter, refCpGList, height, snpPosition, group1);

			epsWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void drawCompactFigure(List<MappedRead> group1, List<MappedRead> group2, int snpPosition,
	                                      String outputFileName) {
		// always put methyl group first
		long methylCountGroup1 = group1.stream()
				.flatMap(r -> r.getCpgList().stream())
				.filter(cpg -> cpg.getMethylStatus() == MethylStatus.C)
				.count();
		long methylCountGroup2 = group2.stream()
				.flatMap(r -> r.getCpgList().stream())
				.filter(cpg -> cpg.getMethylStatus() == MethylStatus.C)
				.count();
		if (methylCountGroup1 < methylCountGroup2) {
			List<MappedRead> tmp = group1;
			group1 = group2;
			group2 = tmp;
		}

		group1 = selectValidReads(group1, snpPosition);
		group2 = selectValidReads(group2, snpPosition);
		List<RefCpG> refCpGList = selectValidSortedRefCpG(group1, group2);
		int imageHeight = (group1.size() + group2.size() + 1) * HEIGHT_INTERVAL, imageWidth = (refCpGList.size() + 3) * CG_RADIUS;
		try {
			EPSWriter epsWriter = new EPSWriter(outputFileName, imageWidth, imageHeight);
			Graphics2D graphWriter = epsWriter.getGraphWriter();
			graphWriter.setBackground(Color.WHITE);
			graphWriter.clearRect(0, 0, imageWidth, imageHeight);
			graphWriter.setPaint(Color.BLACK);
			graphWriter.setFont(new Font("Helvetica", Font.PLAIN, COMMON_FONT_SIZE));

			int height = 0;
			height = drawCompactGroup(graphWriter, refCpGList, height, snpPosition, group1);

			// add line
			graphWriter.setStroke(new BasicStroke(4.0f));
			graphWriter.setPaint(Color.BLUE);
			graphWriter.drawLine(0, height + HEIGHT_INTERVAL / 2, refCpGList.size() * CG_RADIUS,
					height + HEIGHT_INTERVAL / 2);
			graphWriter.setPaint(Color.orange);
			graphWriter.drawLine(refCpGList.size() * CG_RADIUS + CG_RADIUS / 2, 0,
					refCpGList.size() * CG_RADIUS + CG_RADIUS / 2, imageHeight);
			graphWriter.setPaint(Color.BLACK);
			height += HEIGHT_INTERVAL;
			graphWriter.setStroke(new BasicStroke());

			drawCompactGroup(graphWriter, refCpGList, height, snpPosition, group2);

			epsWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<RefCpG> selectValidSortedRefCpG(List<MappedRead> group1, List<MappedRead> group2) {
		Set<RefCpG> validRefCpGSet = new HashSet<>();
		group1.forEach(mappedRead -> mappedRead.getCpgList().forEach(cpg -> validRefCpGSet.add(cpg.getRefCpG())));
		group2.forEach(mappedRead -> mappedRead.getCpgList().forEach(cpg -> validRefCpGSet.add(cpg.getRefCpG())));
		return validRefCpGSet.stream().sorted(RefCpG::compareTo).collect(Collectors.toList());
	}

	private static List<MappedRead> selectValidReads(List<MappedRead> group, int snpPosition) {
		List<MappedRead> validGroup = new ArrayList<>();
		for (MappedRead mappedRead : group) {
			if (snpPosition >= mappedRead.getStart() && snpPosition <= mappedRead.getEnd()) {
				String sequence = mappedRead.getStrand() == '+' ? mappedRead.getSequence() : mappedRead.getComplementarySequence();
				char snp = sequence.charAt(snpPosition - mappedRead.getStart());
				if (snp != '-') {
					validGroup.add(mappedRead);
				}
			}
		}
		return validGroup;
	}

	private static int drawCompactGroup(Graphics2D graphWriter, List<RefCpG> refCpGList, int height, int snpPosition,
	                                    List<MappedRead> group) {
		for (MappedRead mappedRead : group) {
			if (snpPosition >= mappedRead.getStart() && snpPosition <= mappedRead.getEnd()) {
				String sequence = mappedRead.getStrand() == '+' ? mappedRead.getSequence() : mappedRead.getComplementarySequence();
				char snp = sequence.charAt(snpPosition - mappedRead.getStart());
				if (snp != '-') {
					// add cpg
					for (CpG cpG : mappedRead.getCpgList()) {
						if (cpG.getMethylStatus() == MethylStatus.C) {
							graphWriter.fill(
									new Ellipse2D.Double(refCpGList.indexOf(cpG.getRefCpG()) * CG_RADIUS, height,
											CG_RADIUS,
											CG_RADIUS));
						} else if (cpG.getMethylStatus() == MethylStatus.T) {
							graphWriter.draw(
									new Ellipse2D.Double(refCpGList.indexOf(cpG.getRefCpG()) * CG_RADIUS, height,
											CG_RADIUS,
											CG_RADIUS));
						}
					}

					FontMetrics fm = graphWriter.getFontMetrics();
					graphWriter.drawString(String.valueOf(snp),
							(float) ((refCpGList.size() + 1) * CG_RADIUS - fm.charWidth(snp) / 2.0),
							(float) (height + HEIGHT_INTERVAL * 2 / 3.0));
					graphWriter.drawString(String.valueOf(mappedRead.getStrand()),
							(float) ((refCpGList.size() + 2) * CG_RADIUS - fm.charWidth(mappedRead.getStrand()) / 2.0),
							(float) (height + HEIGHT_INTERVAL * 2 / 3.0));
					height += HEIGHT_INTERVAL;
				}
			}
		}
		return height;
	}

	private static void drawFigure(List<MappedRead> group1, List<MappedRead> group2, int minStart, int maxEnd,
	                               int snpPosition, String outputFileName) {
		int imageHeight = (group1.size() + group2.size() + 1) * HEIGHT_INTERVAL, imageWidth = (maxEnd - minStart + 1) * BPWIDTH;
		BufferedImage pngImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graphWriter = pngImage.createGraphics();
		graphWriter.setBackground(Color.WHITE);
		graphWriter.clearRect(0, 0, imageWidth, imageHeight);
		graphWriter.setPaint(Color.BLACK);
		graphWriter.setFont(new Font("Arial", Font.PLAIN, COMMON_FONT_SIZE));

		int height = 0;
		height = drawGroup(graphWriter, minStart, height, snpPosition, group1);

		// add line
		graphWriter.setStroke(new BasicStroke(4.0f));
		graphWriter.setPaint(Color.BLUE);
		graphWriter.drawLine(0, height + HEIGHT_INTERVAL / 2, imageWidth,
				height + HEIGHT_INTERVAL / 2);
		graphWriter.setPaint(Color.BLACK);
		height += HEIGHT_INTERVAL;
		graphWriter.setStroke(new BasicStroke());

		drawGroup(graphWriter, minStart, height, snpPosition, group2);

		File outputPNG = new File(outputFileName);
		try {
			ImageIO.write(pngImage, "png", outputPNG);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static int drawGroup(Graphics2D graphWriter, int minStart, int height, int snpPosition,
	                             List<MappedRead> group) {
		for (MappedRead mappedRead : group) {
			// add line
			graphWriter.drawLine((mappedRead.getStart() - minStart) * BPWIDTH, height + HEIGHT_INTERVAL / 2,
					(mappedRead.getEnd() - minStart + 1) * BPWIDTH, height + HEIGHT_INTERVAL / 2);

			// add cpg
			for (CpG cpG : mappedRead.getCpgList()) {
				if (cpG.getMethylStatus() == MethylStatus.C) {
					graphWriter.fill(
							new Ellipse2D.Double((cpG.getPos() - minStart) * BPWIDTH, height, CG_RADIUS, CG_RADIUS));
				} else if (cpG.getMethylStatus() == MethylStatus.T) {
					graphWriter.draw(
							new Ellipse2D.Double((cpG.getPos() - minStart) * BPWIDTH, height, CG_RADIUS, CG_RADIUS));
				}
			}
			String sequence = mappedRead.getStrand() == '+' ? mappedRead.getSequence() : mappedRead.getComplementarySequence();
			if (snpPosition >= mappedRead.getStart() && snpPosition <= mappedRead.getEnd()) {
				char snp = sequence.charAt(snpPosition - mappedRead.getStart());
				if (snp != '-') {
					graphWriter.drawString(String.valueOf(snp), (snpPosition - minStart) * BPWIDTH,
							height + HEIGHT_INTERVAL * 2 / 3);
				}
			}
			height += HEIGHT_INTERVAL;
		}
		return height;
	}
}
