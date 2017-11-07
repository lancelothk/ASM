ASM-Detector
===
v1.0.0-alpha

# 1.Introduction
ASM-Detector implements a novel computational method to better detect ASMs from whole-genome bisulfite sequencing (WGBS) data. It takes SAM format mapped reads along with FASTA format reference genome as input. In the first step, a serie of parameters are used to select the condidate ASM regions with high quality. Then our graph-based clustering method is applied to detect ASM regions. The output of ASM-Detector includes summary of all detected ASM regions and reads partition results of each ASM region, which can be used in downstream analysis. ASM-Detector also provides several useful tools to perform further analysis of ASM regions, such as visulizing ASM regions in figure.

# 2.Prerequisite
* Java 1.8 or higher

## 2.1 install
Unzip ASM-detector.zip file. Executable files are under bin folder.

To add ASM-detector to $PATH, [check here](http://askubuntu.com/questions/109381/how-to-add-path-of-a-program-to-path-environment-variable)

## 3 Demo
[Download demo data](https://github.com/lancelothk/ASM/releases/download/v1.0.0-alpha/demo.zip)

Demo folder structure:
```
.
├── chr8.fa
└── TRAPPC9_ads_adipose.sam
```
All following command are executed in demo folder. Assuming ASM-Detector required softwares are setup properly(added to $PATH). You can also execute commands through relative path, e.g. ../bin/cpmr, which assume bin folder is in the same level of demo folder.

## 3.1 cpmr: generate candidate ASM regions
```
>cpmr -m TRAPPC9_ads_adipose.sam -r chr8.fa -mcc 4 -mic 5 -mir 10 -p 0.3 --format sam -o TRAPPC9 -pe
Load refChr chr8-0-146364021 with 1309135 RefCpGs. Complete in 2.412000 s
Load 4601 Mapped Reads. Complete in 0.249000 s
```
After execution, demo folder looks like:
```
.
├── chr8.fa
├── TRAPPC9
│   ├── chr8-141108112-141109094.mappedreads
│   ├── chr8-141109226-141110677.mappedreads
│   ├── chr8-141110686-141111080.mappedreads
│   ├── CPMR.bed
│   └── CPMR.report
└── TRAPPC9_ads_adipose.sam
```
## 3.2 asmd: detect ASM regions
```
>asmd -i TRAPPC9 -mic 5 -p 100 -t 8 -o TRAPPC9_result
```
After execution, demo folder looks like:
```
.
├── chr8.fa
├── TRAPPC9
│   ├── chr8-141108112-141109094.mappedreads
│   ├── chr8-141109226-141110677.mappedreads
│   ├── chr8-141110686-141111080.mappedreads
│   ├── CPMR.bed
│   └── CPMR.report
├── TRAPPC9_ads_adipose.sam
└── TRAPPC9_result
    ├── chr8-141108112-141109094.mappedreads.detected
    ├── chr8-141108112-141109094.mappedreads.groups.aligned
    ├── chr8-141109226-141110677.mappedreads.detected
    ├── chr8-141109226-141110677.mappedreads.groups.aligned
    ├── chr8-141110686-141111080.mappedreads.detected
    ├── chr8-141110686-141111080.mappedreads.groups.aligned
    └── detection.summary
```
## 3.3 mfig: generate methylation figure base on aligend reads file.
```
>mfig  -i TRAPPC9_result/chr8-141109226-141110677.mappedreads.groups.aligned -p 141109990 -a T-G -s 22
```
Output is in the same folder of input file with file name extension ".compact.eps".
Example: chr8-141109226-141110677.mappedreads.groups.aligned.compact.eps

# 4.Interface
## 4.1 cpmr:
```
usage: cpmr [options]
    --format <arg>   specify format of input: mappedread, sam
 -m <arg>            MappedRead File (Required)
 -mcc <arg>          Minimum adjacent CpG coverage (Required)
 -mic <arg>          Minimum interval CpG number (Required)
 -mir <arg>          Minimum interval read number (Required)
 -o <arg>            Output Path (Required)
 -p <arg>            Partial methylation threshold (Required)
 -pe                 Pair end mode (Optional)
 -r <arg>            Reference File (Required)
```
## 4.2 asmd:
```
usage: asmd [options]
 -i <arg>     Input intervals folder or interval file name (Required)
 -mic <arg>   Minimum interval CpG number (Required)
 -o <arg>     output folder. Will be create if not exist (Required)
 -p <arg>     Time of random permutation (Required)
 -t <arg>     Thread number to call the program (Required)
 ```
## 4.3 mfig:
```
usage: mfig [options]
 -a <arg>   allele pair. E.g. A-G
 -i <arg>   input grouped read file
 -p <arg>   SNP position
 -s <arg>   font size
```

# 5.Output
## 5.1 cpmr output
### CPMR.bed
Bed file contains every CPMR regions and their read count/cpg count.
### Individual CPMR region read file
Contains reference sequnece and reads with columns(ref,strand,start,end,sequence,ID)
## 5.2 asmd output
### detection.summary
Bed file contains each ASM region detected.
### .detected file
Contains p-value, average methylation information and covered read count(in bracket) of each CpG site in every groups.
### .groups.aligned file
Similar to CPMR region file. However reads in this file are aligned by coordinates and separated into groups for visulization purpose.

---
# Citation
Hu, K., & Li, J. Allele-specific DNA methylation is ubiquitous in human genome and is highly associated with transcriptional regulation (Under Review)

---
# License

Copyright 2017 © Computational Biology lab @ Case Western Reserve University.
See the LICENSE file for license rights and limitations (GPL v3).
