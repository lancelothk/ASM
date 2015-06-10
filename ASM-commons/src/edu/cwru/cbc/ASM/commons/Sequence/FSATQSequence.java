package edu.cwru.cbc.ASM.commons.Sequence;

/**
 * Created by lancelothk on 6/10/15.
 * Represent FASTQ format sequence.
 */
public class FSATQSequence extends Sequence {
	private final String qualSequence;

	public FSATQSequence(String id, String sequence, String qualSequence) {
		super(id, sequence);
		this.qualSequence = qualSequence;
	}
}
