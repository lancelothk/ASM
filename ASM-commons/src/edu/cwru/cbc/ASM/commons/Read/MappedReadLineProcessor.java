package edu.cwru.cbc.ASM.commons.Read;

import edu.cwru.cbc.ASM.commons.CpG.RefCpG;

import java.util.List;

/**
 * Created by ke on 2/21/14.
 * Implemented LineProcessor for reading MappedRead File
 */
public class MappedReadLineProcessor extends MappedReadLineProcessorBase<List<MappedRead>> {
    public MappedReadLineProcessor(List<RefCpG> refCpGList) {
	    super(refCpGList);
    }

	@Override
	public List<MappedRead> getResult() {
        return mappedReadList;
    }
}
