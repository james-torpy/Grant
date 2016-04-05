#!/bin/bash

module load gi/boost/1.53.0
module load gi/samtools/1.2
module load gi/cufflinks/2.2.1

numcores=60

#make directory hierachy
projectname="Grant"
samplename="Grant"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome/annotation files
genomeFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta/GRCh37_iGenome.fa"
annotationFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta/GRCh37_annotation.gtf"

#input/output directories
inType="cufflinks"

inPath="$resultsDir/3.star-cuffdiff_protocol/$samplename.$inType"
outFile="$inPath/all_transcripts"

#scripts/log directory
logDir="$projectDir/scripts/QC/logs"

mkdir -p $logDir

#fetch file names of all transcripts files and put into an array:
i=0
files=( $(ls $inPath/**/transcripts.gtf) )
for file in ${files[@]}; do
	echo -e
	echo The transcripts file used is:
	echo $file
	filesTotal[i]=$file
	let i++;

done;

inFiles=${filesTotal[@]}

echo -e
echo These are the inFiles:
echo $inFiles
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outFile:
echo $outFile
echo -e

#run each transcript.gtf file from cufflinks output for each sample through cuffcompare to merge into one overall transcripts file ready for cuffdiff analysis:
cuffcompare_line="cuffcompare -r $annotationFile -s $genomeFile -o $outFile -CG -R $inFiles"

echo This is the cuffcompare_line:
echo $cuffcompare_line
echo -e

#submit cufflinks job with rame 'CUFFCOMPARE_$samplename' to 16 cluster cores, hold this job until the last cufflinks job is complete:
qsub -N CUFFCOMPARE_$samplename -hold_jid CUFFLINKS_control_proliferative_2 -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $cuffcompare_line
