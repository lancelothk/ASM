package edu.cwru.cbc.ASM.commons.Sequence;

/**
 * Created by lancelothk on 6/10/15.
 * Base class to represent a sequence.
 */
public abstract class Sequence {
	protected final String id;
	protected final String sequence;

	public Sequence(String id, String sequence) {
		this.id = id;
		this.sequence = sequence;
	}

	public String getId() {
		return id;
	}

	public String getSequence() {
		return sequence;
	}

}
