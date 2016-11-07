package edu.cwru.cbc.ASM.commons.io;

import edu.cwru.cbc.ASM.commons.sequence.MappedRead;

import java.io.IOException;

/**
 * Created by kehu on 11/7/16.
 * Handler for processing MappedRead
 */
public interface MappedReadHandler {
	void processMappedRead(MappedRead mappedRead) throws IOException;
}
