#!/bin/bash

module load gi/boost/1.53.0
module load gi/samtools/1.2
module load gi/cufflinks/2.2.1

numcores=60

#make directory hierachy
projectname="Grant"
samplename="Grant"
sample1="control_proliferative_1"
sample2="control_proliferative_2"
sample3="control_proliferative_3"
sample4="case_proliferative_1"
sample5="case_proliferative_2"
sample6="case_proliferative_3"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome file
genomeFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta/GRCh37_iGenome.fa"

echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e

#input/output directories
inType="tophat"
transcripts_ref_inType="cufflinks"
outType="cuffdiff"

inPath="$resultsDir/$samplename.$inType"
transcripts_ref_inPath="$resultsDir/$samplename.$transcripts_ref_inType"
outDir="$resultsDir/$samplename.$outType"

mkdir -p $outDir

echo This is the outDir:
echo $outDir
echo -e

#scripts/log directory
scriptsDir="$projectDir/scripts/QC"
logDir="$scriptsDir/logs"
mkdir -p $logDir

echo This is the logDir:
echo $logDir
echo -e

#fetch transcripts file:
transcripts_inFile=$transcripts_ref_inPath/merged.gtf

echo This is the transcripts_inFile:
echo $transcripts_inFile
echo -e

#input specific samples to be compared
inFiles=( $inPath/$sample1/accepted_hits.bam $inPath/$sample2/accepted_hits.bam $inPath/$sample3/accepted_hits.bam $inPath/$sample4/accepted_hits.bam $inPath/$sample5/accepted_hits.bam $inPath/$sample6/accepted_hits.bam )

echo These are the inFiles:
echo ${inFiles[0]}
echo ${inFiles[1]}
echo ${inFiles[2]}
echo ${inFiles[3]}
echo ${inFiles[4]}
echo ${inFiles[5]}
echo -e

cuffdiff_line="cuffdiff -o $outDir -b $genomeFile -p $numcores -L case_proliferative,control_proliferative -u $transcripts_inFile \
${inFiles[0]},${inFiles[1]},${inFiles[2]} ${inFiles[3]},${inFiles[4]},${inFiles[5]}"

echo $cuffdiff_line

#submit cuffdiff job with rame 'CUFFDIFF_$samplename' to 10 cluster cores, hold this job until CUFFCOMPARE_$uniqueID is complete:
qsub -N CUFFDIFF_$samplename -hold_jid CUFFMERGE_$samplename -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $cuffdiff_line
