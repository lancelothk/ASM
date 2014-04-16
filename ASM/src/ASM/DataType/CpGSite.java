package ASM.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by ke on 2/19/14.
 */
public class CpGSite {
	private long pos;
	private boolean methylated;
    private List<MappedRead> associatedReadList;

	public CpGSite(long pos) {
		this.pos = pos;
		this.methylated = false;
	}

	public CpGSite(long pos, boolean methylated) {
		this.pos = pos;
		this.methylated = methylated;
	}

	public long getPos() {
		return pos;
	}

	public void setPos(long pos) {
		this.pos = pos;
	}

	public boolean isMethylated() {
		return methylated;
	}

	public void setMethylated(boolean methylated) {
		this.methylated = methylated;
	}

    public void addAssociatedRead(MappedRead read) {
        if (associatedReadList == null) {
            associatedReadList = new ArrayList<>();
        }
        associatedReadList.add(read);
    }
}
