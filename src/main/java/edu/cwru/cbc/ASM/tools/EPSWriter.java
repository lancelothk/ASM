package edu.cwru.cbc.ASM.tools;

import net.sf.epsgraphics.ColorMode;
import net.sf.epsgraphics.EpsGraphics;

import java.awt.*;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Created by kehu on 1/6/16.
 */
public class EPSWriter {
	private FileOutputStream epsImage;
	private Graphics2D graphWriter;

	public EPSWriter(String outputFileName, int imageWidth, int imageHeight) throws IOException {
		epsImage = new FileOutputStream(outputFileName + ".eps");
		graphWriter = new EpsGraphics("title", epsImage, 0, 0, imageWidth, imageHeight,
				ColorMode.COLOR_RGB); // eps writer
	}

	public void close() throws IOException {
		((EpsGraphics) graphWriter).close();
		epsImage.close();
		graphWriter.dispose();
	}

	public Graphics2D getGraphWriter() {
		return graphWriter;
	}
}
