options(stringsAsFActors=FALSE)
suppressPackageStartupMessages(library(BitSeq))
suppressPackageStartupMessages(library(data.table))

args=commandArgs(trailingOnly = TRUE)

#get aligned files 
in_file= args[1]
#in_file="/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/Jin/SMART/hg19_cdna_SE/BSF_0137_C64VGACXX_7__Artemether_24h_A6/BSF_0137_C64VGACXX_7__Artemether_24h_A6_bt1_ERCC.sam"
wd= args[2]
sample=unlist(strsplit(basename(in_file),"\\."))[1]
print(sample)

# get reference fasta
#genome=args[2]
ref=args[3]
#ref=paste0("/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/",genome,"/",genome,".fa")
#ref="/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/ERCC92/ERCC92.fa"
print(ref)

# set wd to indikate output directory
setwd(wd)
print(wd)


#Estimate expression

expr=getExpression(in_file,ref,outPrefix=sample,MCMC_chainsN=6,procN=6)
counts=data.table(count=expr$counts)
write.table(counts,paste0(sample,".counts"),sep="\t",row.names=FALSE,quote=FALSE)

print("Done.")
  
