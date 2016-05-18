 #!/bin/bash

module load pethum/bowtie2/prebuilt/2.2.6
module load pethum/tophat/prebuilt/2.1.0
module load gi/samtools/1.2

#number of cores
numcores=8

#genome directories
genomeName="hg38_ercc"

homeDir="/home/jamtor"
genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
bowtie_indexPath="$homeDir/genomes/Homo_sapiens/Ensembl/GRCh37/GenomeFasta/"

#log directory
projectname="Grant"

projectDir="$homeDir/projects/$projectname"
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the bowtie_indexPath:
echo $bowtie_indexPath
echo -e
echo This is the logDir:
echo $logDir
echo -e

#generate the bowtie 2 index files:
bowtie_line="bowtie2-build $genomeFile $genomeName"

echo This is the bowtie_line:
echo $bowtie_line

#submit job with name 'BOWTIE_build_$genomeName' to 10 cluster cores:
qsub -N btbld_$genomeName -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $bowtie_line

#move log files to logDir:
move_line=mv *$genomeDir/$genomeName.o* $logDir

echo This is the move_line:
echo $move_line

#remove uneeded log files:
remove_line=rm *$genomeDir/$genomeName.o*

echo This is the remove_line:
echo $remove_line

#submit above lines to cluster, holding until bowtie build job is done:
qsub -N MOVE_logs -wd $outDir -b y -cwd -j y -R y -pe smp 1 -V $move_line

qsub -N REMOVE_logs -wd $outDir -b y -cwd -j y -R y -pe smp 1 -V $remove_line


