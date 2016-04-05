 #!/bin/bash

module load pethum/bowtie2/prebuilt/2.2.6
module load pethum/tophat/prebuilt/2.1.0
module load gi/samtools/1.2

#number of cores
numcores=10

#genome directories
genomeName="GRCh37_iGenome"
annotationName="GRCh37_annotation"

homeDir="/home/jamtor"
genomeFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta/$genomeName.fa"
annotationFile="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/$annotationName.gtf"
transcriptomeDir="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/tophat.ref"

mkdir -p $transcriptomeDir

#log directory
projectname="Grant"

projectDir="$homeDir/projects/$projectname"
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the transcriptomeDir:
echo $transcriptomeDir
echo -e
echo This is the logDir:
echo $logDir
echo -e

#generate the bowtie 2 index files:
bowtie_line="bowtie2-build -f $genomeFile index"

#generate the tophat transcriptome file:
tophat_trans_line="tophat -G $annotationFile \
--transcriptome-index=$transcriptomeDir \
index"

echo This is the bowtie_line:
echo $bowtie_line
echo -e
echo This is the tophat_trans_line:
echo $tophat_trans_line

#submit job with name 'BOWTIE_build_$genomeName' to 10 cluster cores:
qsub -N BOWTIE_build_$genomeName -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $bowtie_line

#submit job with name 'TOPHAT_trans_$genomeName' to 15 cluster cores:
qsub -N TOPHAT_trans_$genomeName -hold_jid BOWTIE_build_$genomeName -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $tophat_trans_line