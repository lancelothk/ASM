package edu.cwru.cbc.ASM.commons.IO;

import com.google.common.io.LineProcessor;
import edu.cwru.cbc.ASM.commons.Sequence.FASTQSequence;

import java.io.IOException;
import java.util.LinkedHashMap;

/**
 * Created by lancelothk on 6/10/15.
 */
public class FASTQLineProcessor implements LineProcessor<LinkedHashMap<String, FASTQSequence>> {

	@Override public boolean processLine(String line) throws IOException {
		return false;
	}

	@Override public LinkedHashMap<String, FASTQSequence> getResult() {
		return null;
	}
}
