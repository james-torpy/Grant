#load cummeRbund package:
library(cummeRbund)

projectname = "Grant"
genomename = "hg38_ercc"
protocolname = "5.star-cuffdiff_ercc_pairedinfo_protocol"
protocolnumber = "5"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")

inDir= paste0(resultsDir, protocolname, "/Grant.cuffdiff/cummeRbund/")
inFile = paste0(protocolname, "_loaded_data")

scriptDir = paste0(projectDir, "scripts/R")

#set working directory:
setwd(inDir)

#load file containing data:
load(inFile)

writeLines("\n")
paste0("The data loaded is: ", inFile)
writeLines("\n")
paste0("The outDir is: ", inDir)

#create density plot of expression levels in samples and export into .pdf file in current working directory:
filename = paste0(inDir, protocolname, "_densityplot.pdf")
pdf(filename)
csDensity(genes(data))
dev.off()

#set working directory as script directory:
setwd(scriptDir)
