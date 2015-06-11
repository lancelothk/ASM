package edu.cwru.cbc.ASM.commons.Sequence;

/**
 * Created by lancelothk on 6/10/15.
 * Base class to represent a sequence.
 * Since sequence has unique id, id is used to generate hash code.
 */
public abstract class Sequence {
	protected final String id; // should be unique.
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

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;

		Sequence sequence = (Sequence) o;

		return id.equals(sequence.id);

	}

	@Override
	public int hashCode() {
		return id.hashCode();
	}
}
