ASM-Detector
===
v1.0.0-beta

# 1.Introduction
ASM-Detector implements a novel computational method to better detect ASMs from whole-genome bisulfite sequencing (WGBS) data. It takes SAM format mapped reads along with FASTA format reference genome as input. In the first step, a serie of parameters are used to select the condidate ASM regions with high quality. Then our graph-based clustering method is applied to detect ASM regions. The output of ASM-Detector includes summary of all detected ASM regions and reads partition results of each ASM region, which can be used in downstream analysis. ASM-Detector also provides several useful tools to perform further analysis of ASM regions, such as visulizing ASM regions in figure.

# 2.Prerequisite
* Java 1.8 or higher

##2.1 install
Unzip ASM-detector.zip file. Executable files are under bin folder.

To add ASM-detector to $PATH, [check here](http://askubuntu.com/questions/109381/how-to-add-path-of-a-program-to-path-environment-variable)

## 3 Demo


## Citation
Hu, K., & Li, J. Allele-specific DNA methylation is ubiquitous in human genome and is highly associated with transcriptional regulation (Under Review)

## License

Copyright 2017 © Computational Biology lab @ Case Western Reserve University.
See the LICENSE file for license rights and limitations (GPL v3).
