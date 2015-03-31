package edu.cwru.cbc.ASM.commons.DataType;

/**
 * Created by kehu on 11/12/14.
 * Type used to present MethylStatus
 */
public enum MethylStatus {
	C, // methyl
	T, // nonmethyl
	N, // unknown
	E // equal C/T count. Only used in calculating MEC score.
}
