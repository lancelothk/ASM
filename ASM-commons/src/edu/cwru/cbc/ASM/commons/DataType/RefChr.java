package edu.cwru.cbc.ASM.commons.DataType;

/**
 * Created by lancelothk on 5/27/14.
 * Reference chromosome. Contains chr name and reference string.
 */
public class RefChr {
    private String chr;
    private String ref;

    /**
     * @param refString should be upperCase string and without space in it.
     */
    public RefChr(String chr, String refString) {
        this.chr = chr;
        this.ref = refString;
    }

    public String getChr() {
        return chr;
    }

    public String getRefString() {
        return ref;
    }

    public long getStart() {
        return 0;
    }

    public long getEnd() {
        return ref.length() - 1;
    }
}
