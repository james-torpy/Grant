#!/bin/bash

module load gi/samtools/1.2

numcores=10

#make directory hierachy
projectname="Grant"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome directory

genomeDir="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta"
genomeFile="$genomeDir/GRCh37_iGenome.fa"

echo -e
echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e

#input/output types
samplenames=( "Grant" )
inType="trimgalore"
outType="star"

#extension of files to be used:
inExt=".fq.gz"

#scripts/logs directory
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename.$inType"
	outPath="$resultsDir/3.star-cuffdiff_protocol/$samplename.$outType"
	
	echo This is the inPath:
	echo $inPath
	echo -e
	echo This is the outPath:
	echo $outPath
	echo -e

#fetch file names of all projects and put into an array:
	i=0
	files=( $(ls $inPath/**/*$inExt | grep -v unpaired) )
	for file in ${files[@]}; do
		echo The file used is: $file
		echo -e
		filesTotal[i]=$file
		let i++;
	done;

#fetch the inFiles and create an outDir based on their uniqueID:
	j=0
	echo Total files = ${#filesTotal[@]}
	echo -e
	while [ $j -lt ${#filesTotal[@]} ]; do
		inFile1=${filesTotal[$j]}
		inFile2=${filesTotal[$(($j+1))]}
		uniqueID=`basename $inFile1 | sed s/_R1_val_1$inExt//`
		outDir=$outPath/$uniqueID/
			
		mkdir -p $outDir

		echo -e
		echo This is the uniqueID:
		echo $uniqueID
		echo -e
		echo This is the output directory:
		echo $outDir
		echo -e

#align reads of input files with STAR, output into .bam files:
		starJobName="star."$uniqueID
        bamJobName="bam."$uniqueID
        sortJobName="sort."$uniqueID

        indexJobName="index."$uniqueID
        indexStatsJobName="indexstats."$uniqueID
        outSam=$outDir"Aligned.out.sam"
        outBam=$outDir"$uniqueID.bam"
        outSortedBam=$outDir"$uniqueID.sorted.bam"

#in following star_line, the parameters outSAMstrandField and outFilterIntronMotifs make output compatibale with cufflinks:
      	star_line="STAR --runMode alignReads \
       	--readFilesCommand zcat \
     	--genomeDir $genomeDir \
    	--outFilterType BySJout \
      	--outSAMattributes NH HI AS NM MD\
      	--outFilterMultimapNmax 20 \
      	--outFilterMismatchNmax 999 \
      	--outFilterMismatchNoverReadLmax 0.04 \
      	--alignIntronMin 20 \
      	--alignIntronMax 1500000 \
      	--alignMatesGapMax 1500000 \
      	--alignSJoverhangMin 6 \
      	--alignSJDBoverhangMin 1 \
      	--readFilesIn $inFile1 $inFile2 \
      	--outFileNamePrefix $outDir \
      	--runThreadN $numcores \
		--quantMode TranscriptomeSAM \
      	--outFilterMatchNmin 76 \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical"

	echo -e
	echo This is the star_line:
        echo $star_line
        echo -e
        
#remove soft clipping from output .sam file so cufflinks will interpret it correctly:
	remove_sc_line="awk 'BEGIN {OFS="\""\t""\"""} {split(""$""6,C,/[0-9]*/); split(""$""6,L,/[SMDIN]/); \
	if (C[2]=="\""S"\"") {""$""10=substr(""$""10,L[1]+1); ""$""11=substr(""$""11,L[1]+1)}; if (C[length(C)]=="\""S"\"") \
	{L1=length(""$""10)-L[length(L)-1]; ""$""10=substr(""$""10,1,L1); ""$""11=substr(""$""11,1,L1); }; gsub(/[0-9]*S/,"\"""\"",""$""6); print}' \
	$outDir/Aligned.out.sam > $outDir/Aligned.noS.sam"

	echo This is the remove_sc_line:
	echo $remove_sc_line
	echo -e
	
#convert output SAM file from star_line to BAM file using samtools:
	samtools_line="samtools view -bt $genomeFile -o $outDir/$uniqueID.Aligned.out.bam $outDir/Aligned.noS.sam"
	
	echo This is the samtools_line:
	echo $samtools_line
        echo -e

#sort BAM file for cufflinks using samtools:
	BAMsort_line="samtools sort $outDir/$uniqueID.Aligned.out.bam $outDir/$uniqueID.Sorted_Aligned.out"
	
	echo This is the BAMsort_line:
	echo $BAMsort_line
	echo -e

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  	#qsub -N STAR_$uniqueID -hold_jid TRIMGALORE_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $star_line

	#qsub -N REMOVESC_$uniqueID -hold_jid STAR_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp 1 -V $remove_sc_line

	#qsub -N SAMTOOLS_$uniqueID -hold_jid REMOVESC_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $samtools_line

	#qsub -N BAMSORT_$uniqueID -hold_jid SAMTOOLS_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $BAMsort_line

#find out how to do -hold_jid for variable job names!

			j=$(($j+2))

	done;
done;
