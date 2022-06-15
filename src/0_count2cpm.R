setwd("/Users/nanpeng/Documents/projects/ferroptosis/")
########################################################
# step0: load package and function
########################################################
library(edgeR)
library("rtracklayer")
count2cpm = function(count){
  dge<- DGEList(counts=count)
  # 他可以识别count中的非数字？
  keep <- rowSums(dge$counts>0)>0.1*ncol(count)
  #去除的基因考虑进去
  dge <- dge[keep,keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  logcpm <- cpm(dge,log=TRUE) # log2(cpm+1)
  return(logcpm)
}
######################################################
genelist= c("PGD","HMOX1","AKR1C3","FDFT1",
            "FTH1","FANCD2","PTGS2","SQLE","ACSL3","ZEB1",
            "ABCC1","NOX1","CISD1","G6PD","SLC7A11","PEBP1",
            "AIFM2","GLS2","HSBP1","IREB2","GCLC","AKR1C1",
            "NFS1","ACSL4","HSPB1","DPP4","FADS2","PHKG2",
            "RPL8","HMGCR","GSS","TFRC","SLC1A5","CD44","CRYAB",
            "GOT1","STEAP3","ATP5G3","LPCAT3","GCLM","MT1G","ACSF2",
            "EMC2","GPX4","NQO1","ALOX15","CBS","CHAC1","AKR1C2",
            "TP53","KEAP1","NFE2L2","ALOX12","SAT1","ALOX5","CS","ACO1")
######################################################
# step1: read count
######################################################
count.tcga      = read.delim("rawData/tcga/Merge_RNA_seq_Count.txt",row.names = 1,check.names = FALSE)
######################################################
# 1. 
count.gse202203 = data.table::fread("rawData/GSE202203/GSE202203_RawCounts_gene_3207.tsv.gz")
count.gse202203 = data.frame(count.gse202203)
row.names(count.gse202203) = count.gse202203$X
count.gse202203 = count.gse202203[,-1]

genelist[!genelist  %in%  row.names(count.gse202203)]
######################################################
fpkm.gse96058  = data.frame(data.table::fread("rawData/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
row.names(fpkm.gse96058) = fpkm.gse96058$V1
fpkm.gse96058 = fpkm.gse96058[,-1]
genelist[!genelist  %in%  row.names(fpkm.gse96058)]

# fpkm.gse96058.transcript NON log
fpkm.gse96058.transcript  = data.frame(data.table::fread("rawData/GSE96058/GSE96058_transcript_expression_3273_samples_and_136_replicates(1).csv.gz"))
row.names(fpkm.gse96058.transcript) = fpkm.gse96058.transcript$V1
fpkm.gse96058.transcript = fpkm.gse96058.transcript[,-1]
# fpkm.gse96058.transcript 是我们需要的数据
######################################################
geneInfo =  import("rawData/GSE96058/GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz")
geneInfo = geneInfo[geneInfo$type=="exon",]
geneInfo = data.frame(geneInfo)
length = lapply(unique(geneInfo$transcript_id),function(x) {
  sum(geneInfo$width[geneInfo$transcript_id==x])
})
length  = unlist(length)
names(length) = paste(
  unique(geneInfo[,c("transcript_id","geneSymbol")])$transcript_id,
  unique(geneInfo[,c("transcript_id","geneSymbol")])$geneSymbol,sep = "_")
length = data.frame("transcript_id" = gsub("_.*","",names(length)),
                    geneSymbol = gsub(".*_","",names(length)),
                    length = length)
row.names(length) = 1:nrow(length)
saveRDS(length,"result/GSE96058_UCSC_hg38_knownGenes_22sep2014.length.rds")
##############################################################################
row.names(length) = length$transcript_id
fpkm.gse96058.transcript$length = length[row.names(fpkm.gse96058.transcript),"length"]
head(fpkm.gse96058.transcript[,c(1:3,ncol(fpkm.gse96058.transcript))])
temp.count = sweep(fpkm.gse96058.transcript[,1:(ncol(fpkm.gse96058.transcript)-1)],
      STATS=fpkm.gse96058.transcript$length, MARGIN=1,FUN="*")
cpm.gse96058 = temp.count/1000
keep.index = intersect(row.names(length),row.names(cpm.gse96058))
length = length[row.names(length) %in% keep.index,]
cpm.gse96058 = cpm.gse96058[row.names(length),]



cpm.gse96058.gene = cpm.gse96058[!duplicated(length$geneSymbol),]
length.used = length[!duplicated(length$geneSymbol),]
row.names(cpm.gse96058.gene) = length.used$geneSymbol
genelist[!genelist  %in%  row.names(cpm.gse96058.gene)]


# countToFpkm <- function(counts, effLen)
# {
#   N <- sum(counts)
#   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
# }
# fpkm = (count*1e9)/(effLen*N)
# count = (fpkm*effLen*N)/1e9
# cpm   = (count*1e6)/N
# cpm = fpkm/1e3*effLen





######################################################
# step2:count2cpm
######################################################
logcpm.tpm = count2cpm(count.tcga)
logcpm.gse202203 = count2cpm(count.gse202203)
##############################################
# step3:保存
##############################################
write.table(logcpm.tpm,"result/Table/log2cpm_form_count_by_edgeR_tcga.xls",row.names = T,sep = "\t",quote = FALSE)
write.table(logcpm.gse202203,"result/Table/log2cpm_form_count_by_edgeR_gse202203.xls",row.names = T,sep = "\t",quote = FALSE)
write.table(log2(cpm.gse96058.gene+1),"result/Table/log2cpm_form_fpkm_by_gongshi_gse96058.xls",row.names = T,sep = "\t",quote = FALSE)
##############################################





