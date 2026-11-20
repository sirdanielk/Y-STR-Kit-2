@echo off
title Y-STR Kit 1.1

REM     The MIT License (MIT)
REM     Copyright © 2013-2015 Felix Immanuel
REM     http://www.y-str.org
REM
REM     Permission is hereby granted, free of charge, to any person obtaining a copy
REM     of this software and associated documentation files (the Software), to deal
REM     in the Software without restriction, including without limitation the rights
REM     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
REM     copies of the Software, and to permit persons to whom the Software is
REM     furnished to do so, subject to the following conditions: The above copyright
REM     notice and this permission notice shall be included in all copies or substantial
REM     portions of the Software.
REM     THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND.

echo.
echo *** Y-STR Kit 1.1 ***
echo.
echo Tools Used: SAMTools, picard, Cygwin, GATK, Java
echo Project Page: http://www.y-str.org/2015/07/y-str-kit.html
echo Script Developer: Felix Immanuel ^<i@fi.id.au^>
echo.

set "BAM=D:\kls001\bam.hs38.orgfastq_chrY.bam"

echo.
echo Input BAM: %BAM%
echo.

if not exist "%BAM%" (
    echo [ERROR] BAM file not found: %BAM%
    goto ERROR
)

set CYGWIN=nodosfilewarning

echo Pre-execution Cleanup ...

IF EXIST bam_chrY.vcf DEL /Q /F bam_chrY.vcf
IF EXIST bam_chrY.vcf.gz.tbi DEL /Q /F bam_chrY.vcf.gz.tbi
IF EXIST bam_chrY.vcf.idx DEL /Q /F bam_chrY.vcf.idx
IF EXIST bam_complete_sorted.bam DEL /Q /F bam_complete_sorted.bam
IF EXIST bam_complete_sorted.bam.bai DEL /Q /F bam_complete_sorted.bam.bai
IF EXIST bam_sorted.bam DEL /Q /F bam_sorted.bam
IF EXIST bam_sorted.bam.bai DEL /Q /F bam_sorted.bam.bai
IF EXIST bam_sorted_realigned.bai DEL /Q /F bam_sorted_realigned.bai
IF EXIST bam_sorted_realigned.bam DEL /Q /F bam_sorted_realigned.bam
IF EXIST bam_sorted_realigned.bam.bai DEL /Q /F bam_sorted_realigned.bam.bai
IF EXIST bam_wh.bam DEL /Q /F bam_wh.bam
IF EXIST bam_wh_tmp.bam DEL /Q /F bam_wh_tmp.bam
IF EXIST bam_wh_fixed.bam DEL /Q /F bam_wh_fixed.bam
IF EXIST chr DEL /Q /F chr
IF EXIST chrY.bam DEL /Q /F chrY.bam
IF EXIST header DEL /Q /F header
IF EXIST header01 DEL /Q /F header01
IF EXIST header02 DEL /Q /F header02
IF EXIST inchr.bam DEL /Q /F inchr.bam
IF EXIST inchr.sam DEL /Q /F inchr.sam
IF EXIST reads.bam DEL /Q /F reads.bam
IF EXIST tmp.sam DEL /Q /F tmp.sam

echo [OK] Pre-execution Cleanup completed.

echo.
echo Sorting ...
bin\ubin\samtools.exe sort -o bam_complete_sorted.bam "%BAM%"
if errorlevel 1 (
    echo [ERROR] BAM sorting failed.
    goto ERROR
)
if not exist bam_complete_sorted.bam (
    echo [ERROR] bam_complete_sorted.bam was not created.
    goto ERROR
)
echo [OK] Sorting completed.

echo.
echo Indexing the sorted BAM file ...
bin\ubin\samtools.exe index bam_complete_sorted.bam
if errorlevel 1 (
    echo [ERROR] Indexing the sorted BAM file failed.
    goto ERROR
)
if not exist bam_complete_sorted.bam.bai (
    echo [ERROR] BAM index was not created.
    goto ERROR
)
echo [OK] Indexing completed.

echo.
echo Splitting and preparing Chr Y ...
bin\ubin\samtools.exe view -H bam_complete_sorted.bam | bin\ubin\cut -f2 | bin\ubin\grep SN | bin\ubin\grep -E "\b(chr)?Y\b" | bin\ubin\cut -d':' -f2 | bin\ubin\head -1 > chr

if not exist chr (
    echo [ERROR] File chr was not created.
    goto ERROR
)

for /F "tokens=1" %%C in (chr) do (
    bin\ubin\samtools.exe view -b bam_complete_sorted.bam %%C > chrY.bam
)

if not exist chrY.bam (
    echo [ERROR] File chrY.bam was not created.
    goto ERROR
)
for %%F in (chrY.bam) do if %%~zF EQU 0 (
    echo [ERROR] chrY.bam is empty.
    goto ERROR
)
echo [OK] Chromosome Y extraction completed.

echo.
echo Checking and fixing ...
IF EXIST inchr.bam DEL /Q /F inchr.bam
copy chrY.bam inchr.bam > NUL
if errorlevel 1 (
    echo [ERROR] Failed to copy chrY.bam to inchr.bam.
    goto ERROR
)
echo [OK] inchr.bam preparation completed.

echo.
echo Preparing BAM ...
bin\ubin\samtools.exe view inchr.bam | bin\ubin\sed "s/\t/\tchr/2" > tmp.sam
if errorlevel 1 (
    echo [ERROR] Failed to create tmp.sam.
    goto ERROR
)

bin\ubin\cat tmp.sam | bin\ubin\sed "s/\tchrchr/\tchr/" > inchr.sam
if errorlevel 1 (
    echo [ERROR] Failed to create inchr.sam.
    goto ERROR
)

bin\ubin\samtools.exe view -bT ref_hg38\chrY.fa inchr.sam > reads.bam
if errorlevel 1 (
    echo [ERROR] Failed to create reads.bam.
    goto ERROR
)
if not exist reads.bam (
    echo [ERROR] reads.bam was not created.
    goto ERROR
)

bin\ubin\samtools.exe view -H inchr.bam|bin\ubin\grep -v SN > header01
if errorlevel 1 (
    echo [ERROR] Failed to create header01.
    goto ERROR
)

bin\ubin\samtools.exe view -H reads.bam > header02
if errorlevel 1 (
    echo [ERROR] Failed to create header02.
    goto ERROR
)

copy header01+header02 header /Y /B > NUL
if errorlevel 1 (
    echo [ERROR] Failed to merge headers.
    goto ERROR
)

bin\ubin\samtools.exe reheader header reads.bam > bam_wh_tmp.bam
if errorlevel 1 (
    echo [ERROR] Failed to update BAM header.
    goto ERROR
)
echo [OK] BAM preparation completed.

echo Checking BAM header...
bin\ubin\samtools.exe view -H bam_wh_tmp.bam | bin\ubin\grep "^@PG"

echo.
echo Removing conflicts in BAM header ...

bin\ubin\samtools.exe view -H bam_wh_tmp.bam > header
if errorlevel 1 (
    echo [ERROR] Cannot read BAM header.
    goto ERROR
)

bin\ubin\awk.exe -F "\t" "BEGIN { OFS=\"\t\" } /^@PG/ { id=\"\"; for (i=2;i<=NF;i++) if ($i ~ /^ID:/) { id=substr($i,4); break } if (id in seen) next; seen[id]=1 } { print }" header > header01
if errorlevel 1 (
    echo [ERROR] Cannot fix BAM header.
    goto ERROR
)

bin\ubin\samtools.exe reheader -P header01 bam_wh_tmp.bam > bam_wh_fixed.bam
if errorlevel 1 (
    echo [ERROR] Cannot save corrected BAM file.
    goto ERROR
)

bin\ubin\samtools.exe view -H bam_wh_fixed.bam | bin\ubin\grep "^@PG"

echo.
echo Adding read groups and sorting (Picard) ...
bin\jdk17\bin\java.exe -Xmx4g -jar bin\gatk4\gatk-package-4.7.0.0-local.jar AddOrReplaceReadGroups -I bam_wh_fixed.bam -O bam_wh.bam -SORT_ORDER coordinate -RGID rgid -RGLB rglib -RGPL illumina -RGPU rgpu -RGSM sample -VALIDATION_STRINGENCY SILENT

if errorlevel 1 (
    echo [ERROR] Picard AddOrReplaceReadGroups failed.
    goto ERROR
)

if not exist bam_wh.bam (
    echo [ERROR] Picard did not create bam_wh.bam.
    goto ERROR
)

echo [OK] Read groups added successfully.

echo.
echo Sorting ...
bin\ubin\samtools.exe sort -o bam_sorted.bam bam_wh.bam
if errorlevel 1 (
    echo [ERROR] Sorting bam_wh.bam failed.
    goto ERROR
)
if not exist bam_sorted.bam (
    echo [ERROR] bam_sorted.bam was not created.
    goto ERROR
)
echo [OK] Sorting completed.

echo.
echo Indexing the sorted BAM file ...
bin\ubin\samtools.exe index bam_sorted.bam
if errorlevel 1 (
    echo [ERROR] Indexing bam_sorted.bam failed.
    goto ERROR
)
if not exist bam_sorted.bam.bai (
    echo [ERROR] bam_sorted.bam index was not created.
    goto ERROR
)
echo [OK] Indexing completed.

echo.
echo Invoke the variant caller (GATK4 HaplotypeCaller) ...
bin\jdk17\bin\java.exe -Xmx10g -jar bin\gatk4\gatk-package-4.7.0.0-local.jar HaplotypeCaller -R ref_hg38\chrY.fa -I bam_sorted.bam -O bam_chrYc.vcf --ploidy 1 -ERC GVCF --pcr-indel-model NONE --native-pair-hmm-threads 4
if errorlevel 1 (
    echo [ERROR] GATK4 HaplotypeCaller failed.
    goto ERROR
)
if not exist bam_chrYc.vcf (
    echo [ERROR] GATK did not create bam_chrYc.vcf.
    goto ERROR
)
echo [OK] Variant calling completed.

echo.
echo Conversion from GVCF to VCF
bin\jdk17\bin\java.exe -Xmx4g -jar bin\gatk4\gatk-package-4.7.0.0-local.jar GenotypeGVCFs -R ref_hg38\chrY.fa -V bam_chrYc.vcf -O bam_chrY.vcf --sample-ploidy 1 --all-sites
if errorlevel 1 (
    echo [ERROR] Decompression failed.
    goto ERROR
)
if not exist bam_chrY.vcf (
    echo [ERROR] GATK did not create bam_chrY.vcf.
    goto ERROR
)
echo [OK] Decompression completed.


bin\ubin\bgzip.exe bam_chrY.vcf
if errorlevel 1 (
    echo [ERROR] VCF compression failed.
    goto ERROR
)
if not exist bam_chrY.vcf.gz (
    echo [ERROR] bam_chrY.vcf.gz was not created.
    goto ERROR
)
echo [OK] VCF compression completed.

bin\ubin\tabix.exe bam_chrY.vcf.gz
if errorlevel 1 (
    echo [ERROR] VCF indexing with tabix failed.
    goto ERROR
)
if not exist bam_chrY.vcf.gz.tbi (
    echo [ERROR] VCF index was not created.
    goto ERROR
)
echo [OK] VCF indexing completed.

echo.
echo Finding Y-STR Values ...
if not exist out mkdir out
bin\ubin\ystrfinder4.exe ref_hg38\ystr_info.txt bam_chrY.vcf.gz out\Y-STR_Report.html
if errorlevel 1 (
    echo [ERROR] Y-STR Finder failed.
    goto ERROR
)
if not exist out\Y-STR_Report.html (
    echo [ERROR] Y-STR report was not created.
    goto ERROR
)
echo [OK] Y-STR value extraction completed.

REM -- final cleanup
MOVE bam_chrY.vcf.gz out >NUL
if errorlevel 1 (
    echo [ERROR] Failed to move VCF file to the out folder.
    goto ERROR
)
echo [OK] VCF file moved successfully.

IF EXIST bam_chrY.vcf DEL /Q /F bam_chrY.vcf
IF EXIST bam_chrY.vcf.gz.tbi DEL /Q /F bam_chrY.vcf.gz.tbi
IF EXIST bam_chrY.vcf.idx DEL /Q /F bam_chrY.vcf.idx
IF EXIST bam_complete_sorted.bam DEL /Q /F bam_complete_sorted.bam
IF EXIST bam_complete_sorted.bam.bai DEL /Q /F bam_complete_sorted.bam.bai
IF EXIST bam_sorted.bam DEL /Q /F bam_sorted.bam
IF EXIST bam_sorted.bam.bai DEL /Q /F bam_sorted.bam.bai
IF EXIST bam_sorted_realigned.bai DEL /Q /F bam_sorted_realigned.bai
IF EXIST bam_sorted_realigned.bam DEL /Q /F bam_sorted_realigned.bam
IF EXIST bam_sorted_realigned.bam.bai DEL /Q /F bam_sorted_realigned.bam.bai
IF EXIST bam_wh.bam DEL /Q /F bam_wh.bam
IF EXIST bam_wh_tmp.bam DEL /Q /F bam_wh_tmp.bam
IF EXIST bam_wh_fixed.bam DEL /Q /F bam_wh_fixed.bam
IF EXIST chr DEL /Q /F chr
IF EXIST chrY.bam DEL /Q /F chrY.bam
IF EXIST header DEL /Q /F header
IF EXIST header01 DEL /Q /F header01
IF EXIST header02 DEL /Q /F header02
IF EXIST inchr.bam DEL /Q /F inchr.bam
IF EXIST inchr.sam DEL /Q /F inchr.sam
IF EXIST reads.bam DEL /Q /F reads.bam
IF EXIST tmp.sam DEL /Q /F tmp.sam
IF EXIST bam.intervals DEL /Q /F bam.intervals

echo.
out\Y-STR_Report.html

echo All Tasks Completed. Please find results in out subfolder.
echo Also check the logs/info in this window for errors (if any).
goto END
:NOPARAM
echo.
echo  Syntax:
echo     console_bam ^<bam-file^>
echo.
:ERROR
echo.
echo [ERROR] The script failed.
pause
exit /b 1
pause
:END
pause