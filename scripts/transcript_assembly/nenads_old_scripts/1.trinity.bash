#$ -S /bin/bash

module load gi/star/2.3.0e
module unload fabbus/perl/5.14.2
module load fabbus/perl/5.14.2
module load marsmi/bowtie/0.12.8
module load gi/samtools/1.0
module load marsmi/java/1.6.0_37

numcores=24

############## directory hierarchy ##############

#genomeFile="/share/ClusterShare/biodata/contrib/genomeIndices_garvan/pwbc/genomes/Hsapiens/hg19/bowtie2/hg19.fa"
#genomeFile="/share/ClusterShare/biodata/contrib/shinyRseq/hg19_ercc/genome.fa"
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/MRTA"
resultsDir="$projectDir/project_results/"
projectname="ATT_1"

scriptsPath="$projectDir/scripts/transcript_assembly"
logDir=$scriptsPath"/logs"
mkdir -p $logDir


projectname="ATT_1"
inType="trimgalore"

#inPath="$resultsDir/$projectname.$inType"
inPath="/share/ClusterScratch/nenbar/MRTA/project_results/"$projectname.$inType
inPathBit="_R1_001_val_1.fq.gz"

#position of output on clusterScratch
tempFolder="/share/ClusterScratch/nenbar"
tempDir="$tempFolder/MRTA/assembled/$projectname/"
mkdir -p $tempDir

#Get Files
names=($(ls -d $inPath/**/*$inPathBit))

#Set up conditions
for file in ${names[@]};do

	#output    
	shortname=`basename $file | sed 's/-RNA.*_L00//' | sed s/$inPathBit// | sed s/ATT1-//`
	
	fileName=`echo $file | sed s/$inPathBit//`
	echo $fileName
	tempOutPath="$tempDir/$shortname/"
	mkdir -p $tempOutPath

	outDir=$resultsDir/$projectname".trinity"/$shortname
	mkdir -p $outDir

	trinity_version="/home/nenbar/local/lib/test/trinityrnaseq_r20140717/Trinity"
	
	trinity_line="$trinity_version --seqType fq \
		--JM 7G \
		--left $fileName\_R1_001_val_1.fq.gz --right $fileName\_R2_001_val_2.fq.gz \
		--CPU $numcores \
		--SS_lib_type RF \
		--output $tempOutPath \
		--full_cleanup \
		--bflyHeapSpaceMax 6G \
		--normalize_max_read_cov 50 ;"
		#find $tempOutPath -type f -print0 | xargs -0 rm -f"


	trinity_command="trinity_command_$shortname"
	jobName="tri.$shortname"
	echo "#$ -S /bin/bash" >$trinity_command
	echo $trinity_line >>$trinity_command
	chmod 775 $trinity_command
	#echo $trinity_line
	qsub -N $jobName -wd $logDir -j y -R y -pe smp $numcores -V ./$trinity_command

done;

