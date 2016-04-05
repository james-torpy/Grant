#load cummeRbund package:
library(cummeRbund)

projectname = "Grant"
genomename = "hg38_ercc"
protocolname = "5.star-cuffdiff_ercc_pairedinfo_protocol"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
genomeDir = paste0(homeDir, "genomes/", genomename)

inDir = paste0(resultsDir, protocolname, "/Grant.cuffdiff/cummeRbund")
genomeFile = paste0(genomeDir, "/", genomename, ".fa")
annotationFile = paste0(genomeDir, "/gencode.v24.annotation.gtf")

scriptDir = paste0(projectDir, "scripts/R")

system(paste0("mkdir -p", inDir))

#set current working directory:
setwd(inDir)

writeLines("\n")
paste0("The inDir is:")
getwd()
writeLines("\n")
paste0("The genomeFile is: ", genomeFile)
writeLines("\n")
paste0("The annotationFile is: ", annotationFile)

#load cuffdiff data with genome fasta file and annotation file:
data=readCufflinks(genome = genomeFile, gtfFile = annotationFile, rebuild=TRUE)

#save data as image:
save.image(paste0(protocolname, "_loaded_data.Rdata"))

#set working directory as script directory:
setwd(scriptDir)