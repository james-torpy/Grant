#!/bin/bash

module load pethum/bowtie2/prebuilt/2.2.6
module load pethum/tophat/prebuilt/2.1.0
module load gi/samtools/1.2

numcores=32

#make directory hierachy
projectname="neura"

homeDir="/share/Temp/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome/reference transcriptome directory
genomeName="GRCh37_iGenome"
annotationName="GRCh37_annotation.gtf"

genomeDir="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta"
genomeFile="$genomeDir/$genomeName"
annotationFile="$genomeDir/$annotationName"

echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile

#input/output types
samplenames=( "samples" )
inType="trimgalore"
outType="tophat"

#extension of files to be used:
inExt=".fq"

#scripts/logs directory
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename.$inType"
	outPath="$resultsDir/$samplename.$outType"
	
	echo This is the inPath:
	echo $inPath
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
		echo This is the outDir:
		echo $outDir
		echo -e

#align reads of input files with tophat, output into .bam files:
		tophat_line="tophat -p $numcores -G $annotationFile -o $outDir $genomeFile $inFile1 $inFile2"

		echo This is the tophat_line:
		echo $tophat_line
		echo -e

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
  		qsub -N TOPHAT_$uniqueID -hold_jid TRIMGALORE_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $tophat_line

		j=$(($j+2))

	done;
done;
