# Y-STR Kit 2 Bundle
- This is based on original work of [Y-STR Kit](https://github.com/savolab/Y-STR_Kit) author, to which credits go on
- This bundle contains original application with applied improvements (ystrfinder)
- Plus another application adjusted to work with GATK4 (ystrfinder4)

## 1. Purpose of this project
- Main goal was to adapt code to work with newer GATK4 HaplotypeCaller, since GATKv3 genotyping tools are considered obsolete

## 2. What's new
- Introduced two more reference ref genomes into workflow (hg38 and t2t), apart from existing hg37
- Added support for GATK4 HaplotypeCaller
- Modified regex patterns for motifs and filters to make output match my existing data
- Improved handling of pesky bam headers that exist often in ancient genome samples
- Improved html viewer with added colours - red for deletions and green for insertions, no calls or unreliable calls stays grey (GATK4 pipeline has no green because of different format of gVCF)

## 3. Before first use
- set input files in .bat files
- place required tools in appropriate subfolders
- build ref genomes if not provided

## 4. The outcome
- There is no reliable single solution STR genotyping tool. The cause is that output may rely on aligner work and version of ref genome used.
  Also many regions of chromosome Y are hard to map. That way, adding another pipeline does not resolve all problems. Instead, after some testing the idea would be like to:
  - use both new reference genomes and new gatk4 pipeline to get four results from four different runs and compare them statistically, choosing the most frequent value
  - check the html viewer to find which of reads are clean or which reads can be corrected based on visible DNA string below STR table

## 5. To do
- Fine tune and correct some motifs and regex
- Get more use of information inside gVCF format
- Finally leave slow java based apps to make new python based workflow
- Find more reliable automated STR genotyping solution in python

## 6. 
- This project is based on the original idea and work by original Author [Y-STR Kit](https://github.com/savolab/Y-STR_Kit). The software is provided as is and the licence is the same as the original. 

## 7. AI Disclaimer
- AI Chat GPT was used to resolve errors appeared after switching to GATK4, to improve handling of ancient genomes, to improve html viewer.
