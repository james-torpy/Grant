#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=15

#genome directories
genome_fileName="GRCh37_iGenome.fa"
annotation_fileName="GRCh37_annotation.gtf"

homeDir="/home/jamtor"
genomeDir="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta"
genomeFile="$genomeDir/$genome_fileName"
annotationFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/$annotation_fileName"
outDir="$genomeDir/../star.ref/"

mkdir -p $outDir

scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir
echo $outDir
echo -e

#generate the star reference files:
star_ref_line="$homeDir/local/bin/STAR \
	--runMode genomeGenerate \
	--sjdbGTFfile $annotationFile --sjdbOverhang 99 --genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --runThreadN $numcores \
	--outFileNamePrefix $outDir"

echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genome -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $star_ref_line