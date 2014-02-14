package ASM;

import java.util.Comparator;

public class ReadComparator implements Comparator<Read> {

    @Override
    public int compare(Read r1, Read r2) {
        // rank smaller start, longer end in front
        if (r1.getStart() < r2.getStart()) {
            return -1;
        } else if (r1.getStart() > r2.getStart()) {
            return 1;
        } else if (r1.getEnd() < r2.getEnd()) {
            return 1;
        } else if (r1.getEnd() > r2.getEnd()) {
            return -1;
        } else {
            return 0;
        }
    }
}
