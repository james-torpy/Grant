#load cummeRbund package:
library(cummeRbund)

projectname = "Grant"
genomename = "hg38_ercc"
protocolname = "3.star-cuffdiff_protocol"
protocolnumber = "3"

homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
genomeDir = paste0(homeDir, "genomes/", genomename)

inDir= paste0(resultsDir, protocolname, "/Grant.cuffdiff/")
genomeFile=paste0(genomeDir, "/", genomename, ".fa")
annotationFile=paste0(genomeDir, "/gencode.v24.annotation.gtf")
outPrefix=paste0(inDir, "/Grant_protocol_", protocolnumber, "_")

#set current working directory:
setwd(inDir)
getwd()

#load cuffdiff data with genome fasta file and annotation file:
Grant_data=readCufflinks(genome = genomeFile, gtfFile = annotationFile, rebuild=TRUE)

#create density plot of expression levels in samples and export into .png file in current working directory
csDensity(genes(Grant_data))

#open .png file with desired filename
filename = paste0(outPrefix,"densityplot.png")
filename
png(filename = filename)
    
csDensity(genes(Grant_data))
dev.off()

#create scatterplot to compare expression of each gene in two conditions and export into .png file in current working directory:
csScatter(genes(Grant_data), 'control_proliferative', 'case_proliferative')

filename = paste0(outPrefix,"scatterplot.png")
filename
png(filename = filename)
csScatter(genes(Grant_data), 'control_proliferative', 'case_proliferative')
dev.off()

#create volcano plot to inspect differentially expressed genes and export into .png file in current working directory:
csVolcano(genes(Grant_data), 'control_proliferative', 'case_proliferative')

filename = paste0(outPrefix,"volcanoplot.png")
filename
png(filename = filename)
csVolcano(genes(Grant_data), 'control_proliferative', 'case_proliferative')
dev.off()

#specify gene of interest:
gene_of_interest = getGene(Grant_data, 'GAPDH')

#plot expression levels for gene of interest on bar plot and export into .png file in current working directory:
expressionBarplot(gene_of_interest, logMode=T)

filename = paste0(outPrefix,"GAPDH_barplot.png")
filename
png(filename = filename)
expressionBarplot(gene_of_interest, logMode=T)
dev.off()

#plot expression levels for isofroms for gene of interest on bar plot and export into .png file in current working directory:
expressionBarplot(isoforms(gene_of_interest), logMode=T)

filename = paste0(outPrefix,"GAPDH_isoform_barplot.png")
filename
png(filename = filename)
expressionBarplot(isoforms(gene_of_interest), logMode=T)
dev.off()

save.image(paste0("/Users/jamestorpy/clusterHome/projects/Grant/results/", protocolname, "/Grant.cuffdiff/Grant_protocol_", protocolnumber, ".Rdata"))
