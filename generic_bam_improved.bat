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

set "BAM=D:\kls001\bam.hs19.orgfastq_chrY.bam"

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
bin\ubin\samtools.exe sort --no-PG -o bam_complete_sorted.bam "%BAM%"
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
    echo [ERROR] Indexing the sorted BAM failed.
    goto ERROR
)
if not exist bam_complete_sorted.bam.bai (
    echo [ERROR] BAM index was not created.
    goto ERROR
)
echo [OK] Indexing completed.

echo.
echo Splitting and preparing Chr Y ...
bin\ubin\samtools.exe view -H bam_complete_sorted.bam|bin\ubin\cut -f2|bin\ubin\grep SN|bin\ubin\grep Y|bin\ubin\cut -d':' -f2|bin\ubin\head -1 > chr

if not exist chr (
    echo [ERROR] chr file was not created.
    goto ERROR
)

for /F "tokens=1" %%C in (chr) do (
    bin\ubin\samtools.exe view -b bam_complete_sorted.bam %%C > chrY.bam
)

if not exist chrY.bam (
    echo [ERROR] chrY.bam was not created.
    goto ERROR
)
for %%F in (chrY.bam) do if %%~zF EQU 0 (
    echo [ERROR] chrY.bam is empty.
    goto ERROR
)
echo [OK] Y chromosome extraction completed.

echo.
echo Checking and fixing ...
IF EXIST inchr.bam DEL /Q /F inchr.bam
copy chrY.bam inchr.bam > NUL
if errorlevel 1 (
    echo [ERROR] Copying chrY.bam to inchr.bam failed.
    goto ERROR
)
echo [OK] inchr.bam preparation completed.

echo.
echo Preparing BAM ...
bin\ubin\samtools.exe view inchr.bam | bin\ubin\sed "s/\t/\tchr/2" > tmp.sam
if errorlevel 1 (
    echo [ERROR] Creating tmp.sam failed.
    goto ERROR
)

bin\ubin\cat tmp.sam | bin\ubin\sed "s/\tchrchr/\tchr/" > inchr.sam
if errorlevel 1 (
    echo [ERROR] Creating inchr.sam failed.
    goto ERROR
)

bin\ubin\samtools.exe view -bT ref\chrY.fa inchr.sam > reads.bam
if errorlevel 1 (
    echo [ERROR] Creating reads.bam failed.
    goto ERROR
)
if not exist reads.bam (
    echo [ERROR] reads.bam was not created.
    goto ERROR
)

bin\ubin\samtools.exe view -H inchr.bam|bin\ubin\grep -v SN > header01
if errorlevel 1 (
    echo [ERROR] Creating header01 failed.
    goto ERROR
)

bin\ubin\samtools.exe view -H reads.bam > header02
if errorlevel 1 (
    echo [ERROR] Creating header02 failed.
    goto ERROR
)

copy header01+header02 header /Y /B > NUL
if errorlevel 1 (
    echo [ERROR] Merging the headers failed.
    goto ERROR
)

bin\ubin\samtools.exe reheader header reads.bam > bam_wh_tmp.bam
if errorlevel 1 (
    echo [ERROR] Updating the BAM header failed.
    goto ERROR
)
echo [OK] BAM preparation completed.

echo Checking the BAM header...
bin\ubin\samtools.exe view -H bam_wh_tmp.bam | bin\ubin\grep "^@PG"

echo.
echo Removing conflicts in the BAM header ...

bin\ubin\samtools.exe view -H bam_wh_tmp.bam > header
if errorlevel 1 (
    echo [ERROR] Unable to read the BAM header.
    goto ERROR
)

bin\ubin\awk.exe -F "\t" "BEGIN { OFS=\"\t\" } /^@PG/ { id=\"\"; for (i=2;i<=NF;i++) if ($i ~ /^ID:/) { id=substr($i,4); break } if (id in seen) next; seen[id]=1 } { print }" header > header01
if errorlevel 1 (
    echo [ERROR] Unable to fix the BAM header.
    goto ERROR
)

bin\ubin\samtools.exe reheader -P header01 bam_wh_tmp.bam > bam_wh_fixed.bam
if errorlevel 1 (
    echo [ERROR] Unable to save the corrected BAM file.
    goto ERROR
)

bin\ubin\samtools.exe view -H bam_wh_fixed.bam | bin\ubin\grep "^@PG"

echo.
echo Adding read groups ...

bin\jre\bin\java.exe -Xmx2g -jar bin\picard\picard.jar AddOrReplaceReadGroups -INPUT bam_wh_fixed.bam -OUTPUT bam_wh.bam -SORT_ORDER coordinate -RGID rgid -RGLB rglib -RGPL illumina -RGPU rgpu -RGSM sample -VALIDATION_STRINGENCY SILENT

if errorlevel 1 (
    echo [ERROR] Picard AddOrReplaceReadGroups failed.
    goto ERROR
)

if not exist bam_wh.bam (
    echo [ERROR] Picard did not create bam_wh.bam.
    goto ERROR
)

echo [OK] Adding read groups completed.

echo.
echo Sorting ...
bin\ubin\samtools.exe sort --no-PG -o bam_sorted.bam bam_wh.bam
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
    echo [ERROR] The bam_sorted.bam index was not created.
    goto ERROR
)
echo [OK] Indexing completed.

echo.
echo Realignment of the sorted and indexed BAM file ...
bin\jre\bin\java.exe -Xmx2g -jar bin\gatk\GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref\chrY.fa -I bam_sorted.bam -o bam.intervals
if errorlevel 1 (
    echo [ERROR] GATK RealignerTargetCreator failed.
    goto ERROR
)
if not exist bam.intervals (
    echo [ERROR] GATK did not create bam.intervals.
    goto ERROR
)
echo [OK] RealignerTargetCreator completed.

bin\jre\bin\java.exe -Xmx2g -jar bin\gatk\GenomeAnalysisTK.jar -T IndelRealigner -R ref\chrY.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam
if errorlevel 1 (
    echo [ERROR] GATK IndelRealigner failed.
    goto ERROR
)
if not exist bam_sorted_realigned.bam (
    echo [ERROR] GATK did not create bam_sorted_realigned.bam.
    goto ERROR
)
echo [OK] Realignment completed.

echo.
echo Indexing the realigned BAM file ...
bin\ubin\samtools.exe index bam_sorted_realigned.bam
if errorlevel 1 (
    echo [ERROR] Indexing bam_sorted_realigned.bam failed.
    goto ERROR
)
if not exist bam_sorted_realigned.bam.bai (
    echo [ERROR] The realigned BAM index was not created.
    goto ERROR
)
echo [OK] Indexing the realigned BAM completed.

echo.
echo Invoke the variant caller ...
bin\jre\bin\java.exe -Xmx2g -jar bin\gatk\GenomeAnalysisTK.jar -l INFO -R ref\chrY.fa -T UnifiedGenotyper -glm BOTH -I bam_sorted_realigned.bam -rf BadCigar -nct 2 -o bam_chrY.vcf --output_mode EMIT_ALL_CONFIDENT_SITES
if errorlevel 1 (
    echo [ERROR] GATK UnifiedGenotyper failed.
    goto ERROR
)
if not exist bam_chrY.vcf (
    echo [ERROR] GATK did not create bam_chrY.vcf.
    goto ERROR
)
echo [OK] Variant calling completed.

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
    echo [ERROR] Indexing the VCF with tabix failed.
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
bin\ubin\ystrfinder.exe ref\ystr_info.txt bam_chrY.vcf.gz out\Y-STR_Report.html
if errorlevel 1 (
    echo [ERROR] Y-STR Finder failed.
    goto ERROR
)
if not exist out\Y-STR_Report.html (
    echo [ERROR] Y-STR report was not created.
    goto ERROR
)
echo [OK] Reading Y-STR values completed.

REM -- final cleanup
MOVE bam_chrY.vcf.gz out >NUL
if errorlevel 1 (
    echo [ERROR] Moving the VCF file to the out folder failed.
    goto ERROR
)
echo [OK] Moving the VCF file completed.

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
echo.
:ERROR
echo.
echo [ERROR] The script failed.
pause
exit /b 1
pause
:END
pause
