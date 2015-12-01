package edu.cwru.cbc.ASM.commons.methylation;

import net.openhft.koloboke.collect.map.hash.HashIntObjMap;
import net.openhft.koloboke.collect.map.hash.HashIntObjMaps;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by kehu on 11/13/14.
 * MethylationUtils for methods shared by modules under ASM project.
 */
public class MethylationUtils {
	public static final int REFERENCE_INIT_POS = 0;

	/**
	 * Extract reference CpG from reference string.
	 *
	 * @param initPos   0-based. Diff to UCSCGB, which use 1-based start and end.
	 * @param reference should be upperCase string and without space in it.
	 * @return returned List is sorted by position.
	 */
	public static List<RefCpG> extractCpGSite(String reference, int initPos) {
		// initialize refCpGList with 1/20 of reference size, since the probability of CG occurrence should be 1/16.
		// Actual probability is higher than 1/16 from current observation. In CPG dense region, it should be much higher.
		List<RefCpG> reFCpGList = new ArrayList<>(reference.length() / 20);
		for (int i = 0; i < reference.length() - 1; i++) {
			if (reference.charAt(i) == 'C' && reference.charAt(i + 1) == 'G') {
				reFCpGList.add(new RefCpG(i + initPos));
			}
		}
		return reFCpGList;
	}

	public static HashIntObjMap<RefCpG> initialzeRefCpGMap(RefChr refChr) {
		List<RefCpG> refCpGList = extractCpGSite(refChr.getRefString(), MethylationUtils.REFERENCE_INIT_POS);
		HashIntObjMap<RefCpG> refMap = HashIntObjMaps.newMutableMap();
		for (RefCpG refCpG : refCpGList) {
			refMap.put(refCpG.getPos(), refCpG);
		}
		return refMap;
	}

	public static Map<String, HashIntObjMap<RefCpG>> initializeGenomeRefCpGMap(Map<String, RefChr> genomeReferenceMap) {
		Map<String, HashIntObjMap<RefCpG>> genomeRefCpGMap = new HashMap<>();
		for (Map.Entry<String, RefChr> entry : genomeReferenceMap.entrySet()) {
			genomeRefCpGMap.put(entry.getKey(), initialzeRefCpGMap(entry.getValue()));
		}
		return genomeRefCpGMap;
	}
}
