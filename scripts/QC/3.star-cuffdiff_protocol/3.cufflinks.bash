#!/bin/bash

module load gi/boost/1.53.0
module load gi/samtools/1.2
module load gi/cufflinks/2.2.1

numcores=10

#make directory hierachy
projectname="Grant"
samplename="Grant"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output directories
inType="star"
outType="cufflinks"

inPath="$resultsDir/$samplename.$inType"
outPath="$resultsDir/3.star-cuffdiff_protocol/$samplename.$outType"

#scripts/log directory
scriptsDir="$projectDir/scripts/QC"
logDir="$scriptsDir/logs"

mkdir -p $logDir

#fetch directory names for files:
i=0
directory_names=`ls $inPath`
for directory in ${directory_names[@]};do
        echo The directory used is: $directory;
        echo -e
        directoriesTotal[i]=$directory
        let i++
done;

echo The total number of directories is: ${#directoriesTotal[@]}
echo -e

#set up conditions to perform analysis on bam files from STAR alignment output

#set up directories specific to each file being analysed:
j=0

while [ $j -lt ${#directoriesTotal[@]} ]; do
        
	uniqueID=`basename ${directoriesTotal[$j]}`
        unique_inDir=$inPath/${directoriesTotal[$j]}
        inFile=$unique_inDir/$uniqueID.Sorted_Aligned.out.bam
		uniqueID=`basename $unique_inDir`
        outDir=$outPath/${directoriesTotal[$j]}
 
	mkdir -p $outDir

        echo -e
        echo This is the inFile:
        echo $inFile
        echo -e
	echo This is the outDir:
	echo $outDir
	echo -e
	echo -e This is the uniqueID:
	echo $uniqueID
	echo -e

#perform analysic on bam files using the reference genome and annotation file:
	cufflinks_line="cufflinks -p $numcores -o $outDir $inFile"

	echo This is the cufflinks_line:
	echo $cufflinks_line
	echo -e

#submit cufflinks job with rame 'CUFFLINKS_$uniqueID' to 10 cluster cores, hold this job until STAR_$uniqueID is complete:
	qsub -N CUFFLINKS_$uniqueID -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $cufflinks_line

	j=$(($j+1))

done;
