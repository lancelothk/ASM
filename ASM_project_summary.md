ASM project summary
=======

1. Problem and Goal:
-----
##### Data: 
Bisulfite sequencing data, CG/CT as methylation infomation.   

##### Goals: 
+ Detect ASM region genome-widely
    + select condidate regions
    + check is the region is ASM region -- the region with two clusters of reads which contain different methylation pattern
    + assemble two distinct methylation patterns which come from different alleles.

##### Formal notation:
* Input: A set of fragments R in which each fragment is consist of three values {C,T,N}. The fragments are pre-aligned to a reference fragment(a long fragment with continuous sites). Each of the character in the input fragment covers a site in the reference fragment. 
* Output: 
    * A partation P(R1, R2) in which R1 and R2 are set of input fragments. Two assembled methylation patterns, M1 from R1, M2 from R2. M1 and M2 should be different. 
    * If R1/R2 is empty or R1 and R2 have same/simillar methylation pattern or |R1| >> |R2| or |R2| >> |R1|, this region is not a ASM region. Otherwise, this region is a detected ASM region.


##### Data preprocessing:
  1. Only use reads covers at least two CpG sites. 
  2. Only CpG with coverage > MIN_CONT_COVERAGE is included in candidate interval.
  3. Only intervals with reads > MIN_INTERVAL_READ and CpG number > MIN_INTERVAL_CPG is included in following analysis.


2. Cons of existing methods:
-----

3. Methods:
-----

4. Result evaluation:
-----
