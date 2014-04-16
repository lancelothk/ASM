package ASM.DataType;

import java.util.Comparator;

/**
 * Created by ke on 4/9/14.
 */
public class MappedReadComparaterByCpG implements Comparator<MappedRead> {
    @Override
    public int compare(MappedRead o1, MappedRead o2) {
        if (o1.getFirstCpG().getPos() == o2.getFirstCpG().getPos()) {
            return (int) (o1.getLastCpG().getPos() - o2.getLastCpG().getPos());
        } else {
            return (int) (o1.getFirstCpG().getPos() - o2.getFirstCpG().getPos());
        }
    }
}
