project.init()
inputFolder = "/fhgfs/groups/lab_bock/shared/COREseq/results_pipeline"

pipeDirs = list.dirs(inputFolder, recursive=FALSE)

results=list()
dir = pipeDirs[[1]];
for (dir in pipeDirs) {
	message(dir);
	statFiles = list.files(dir, pattern="_stats.txt")
	statFiles2 = list.files(dir, pattern="stats_")
	statFiles = c(statFiles, statFiles2)
	for (statFile in statFiles) {
		message(statFile);
		pipeline = gsub("_stats.txt", "", statFile)
		pipeline = gsub("stats_", "", pipeline)
		statPath = paste0(dir, "/", statFile);
		a = fread(statPath)
		setnames(a, c("key", "value"))
		a[,key:=gsub(" ", "_", key)] # Change spaces to underscores
		setkey(a, "key")
		a[,sampleName:=basename(dir)]
		a[,pipeline:=pipeline]
		sampleName = basename(dir)
		results[[sampleName]] = a;
	}
}
resultsDT = do.call(rbind, results)
# Select latest for identical statistics
resultsDT = resultsDT[,list(value=value[length(value)]), by=c("key", "sampleName", "pipeline"), roll=+Inf]


# clean out garbage columns
garbageCols = c("K1_unmethylated_CHG",
"K1_unmethylated_CHG_count",
"K1_unmethylated_CHG_pct",
"K1_unmethylated_CHH",
"K1_unmethylated_CHH_count",
"K1_unmethylated_CHH_pct",
"K1_unmethylated_CpG",
"K1_unmethylated_CpG_count",
"K1_unmethylated_CpG_pct",
"K3_methylated_C_count" ,
"K3_methylated_CHG",
"K3_methylated_CHG_count",
"K3_methylated_CHG_pct",
"K3_methylated_CHH",
"K3_methylated_CHH_count",
"K3_methylated_CHH_pct",
"K3_methylated_CpG",
"K3_methylated_CpG_count",
"K3_methylated_CpG_pct"  )

resultsDT = resultsDT[! key %in% garbageCols,]

library(reshape2) #no longer necessary after data.table 1.9.5??
resultsTable = dcast(resultsDT, formula = "... ~ key")
write.tsv(resultsTable, paste0(inputFolder, "/qc_results_complete.tsv"))

# Make table prettier
resultsTable = as.data.table(resultsTable)
class(resultsTable)
allCols = colnames(resultsTable)
allCols
resultsTable
resultsTable[1:10,]


resultsTable[, trim_rate := (Bam_reads - Trimmed_size)/Bam_reads]
resultsTable[, alignment_rate := (Trimmed_size - Aligned_reads)/Trimmed_size]
resultsTable[, dupe_rate := (Aligned_reads - Deduplicated)/Aligned_reads]
resultsTable[, filt_rate := (Deduplicated - Filtered)/Deduplicated]
resultsTable

resultsTable[pipeline=="WGBS",  c("sampleName", "pipeline", "Bam_reads", "Trimmed_size", "trim_rate", "Aligned_reads", "alignment_rate", "Deduplicated", "dupe_rate", "Filtered", "filt_rate"), with=F]

setcolorder(resultsTable, c("sampleName", "pipeline", "Bam_reads", "Trimmed_size", "trim_rate", "Aligned_reads", "alignment_rate", "Deduplicated", "dupe_rate", "Filtered", "filt_rate", ))

