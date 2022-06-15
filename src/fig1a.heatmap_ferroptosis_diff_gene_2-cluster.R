######################################################
# step0. load package and function
######################################################
library(ComplexHeatmap)
library(TCGAbiolinks)
library(edgeR)
library(limma)
scale_rows=function(x){
  m=apply(x,1,mean,na.rm=T)
  s=apply(x,1,sd,na.rm=T)
  return((x-m)/s)
}
scale_mat=function(mat,scale){
  if(!(scale%in%c("none","row","columa"))){
    stop("scale argument should take value:'none','row'or'column'")
  }
  mat = switch(scale,none=mat,row=scale_rows(mat),column=
                 t(scale_rows(t(mat))))
  return(mat)
}
######################################################
setwd("/Users/nanpeng/Documents/projects/ferroptosis/")
######################################################
count = read.delim("rawdata/BRAC/Merge_RNA_seq_Count.txt",row.names = 1,check.names = FALSE)
# strsplit("TCGA-H6-8124-11",split = "-")[[1]][4]
matadata = data.frame(
  "sampleID" = colnames(count),
  "group" = ifelse(as.numeric(gsub(".*-","",colnames(count))) > 10,"normal","tumor")
)

table(matadata$group)
######################################################
# step2: DEGs of 铁死亡-related genes
######################################################
# step2.1: 给定铁死亡-related genes
######################################################
genelist= c("ENSG00000164120","ENSG00000100292","ENSG00000196139","ENSG00000079459",
            "ENSG00000167996","ENSG00000144554","ENSG00000073756","ENSG00000104549",
            "ENSG00000123983","ENSG00000148516","ENSG00000103222","ENSG00000120658",
            "ENSG00000122873","ENSG00000160211","ENSG00000151012","ENSG00000089220",
            "ENSG00000042286","ENSG00000135423","ENSG00000230989","ENSG00000136381",
            "ENSG00000001084","ENSG00000187134","ENSG00000244005","ENSG00000068366",
            "ENSG00000106211","ENSG00000197635","ENSG00000134824","ENSG00000156873",
            "ENSG00000161016","ENSG00000113161","ENSG00000100983","ENSG00000072274",
            "ENSG00000105281","ENSG00000026508","ENSG00000109846","ENSG00000120053",
            "ENSG00000115107","ENSG00000154518","ENSG00000111684","ENSG00000023909",
            "ENSG00000125144","ENSG00000167107","ENSG00000104412","ENSG00000167468",
            "ENSG00000181019","ENSG00000161905","ENSG00000160200","ENSG00000128965",
            "ENSG00000151632","ENSG00000141510","ENSG00000079999","ENSG00000116044",
            "ENSG00000108839","ENSG00000212766","ENSG00000012779","ENSG00000062485",
            "ENSG00000204178"
)
######################################################
# step2.2: 计算基因表达量，并提取铁死亡-related genes表达量
######################################################
dge<- DGEList(counts=count,group = matadata$group)
# 他可以识别count中的非数字？
keep <- rowSums(dge$counts>0)>0.1*ncol(count)
#去除的基因考虑进去
dge <- dge[keep,keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logcpm <- cpm(dge,log=TRUE) # log2(cpm+1)
###############################################################
#（可选）求pca，去除偏离值
# install.packages("factoextra")
library(factoextra)
pca_analysis = function(fpkm,metadata){
  varValue = apply(fpkm,1,var)
  varValue = sort(varValue,decreasing = T)
  expr_mat = fpkm[head(names(varValue),5000),]
  print(tail(head(varValue,10000)))
  if(max(expr_mat)>50){
    data.pca = prcomp(t(log2(expr_mat+1)), center = TRUE,scale = TRUE)
  }else{
    data.pca = prcomp(t(expr_mat), center = TRUE,scale = TRUE)
  }
  p1 = fviz_pca_ind(data.pca,addEllipses=TRUE, ellipse.level=0.99,label = "none")
  # palette = c("#dfc27d","#00AFBB"),
  #habillage=metadata[gsub("\\.","-",row.names(data.pca$x)),"TME2.subtype"])
  result = list("pcaMat" = data.pca,"pcaPlot" = p1)
  return(result)
}
plotpca = pca_analysis(logcpm,annotation.col)
plotpca$pcaPlot
pcaMat = data.frame(plotpca$pcaMat$x)
row.names(pcaMat)[!((pcaMat$PC2>=50 &pcaMat$PC1>=40 ) | (pcaMat$PC1>60 &pcaMat$PC2>=25 ))]
################################################################
logcpm= logcpm[,row.names(pcaMat)[!((pcaMat$PC2>=30 &pcaMat$PC1<=-40 ) | (pcaMat$PC2>=45 &pcaMat$PC1<=0 ))]]
################################################################
# exprMat是铁死亡-related genes表达量
exprMat=logcpm[rownames(logcpm) %in% genelist,]
################################################################
# step2.3: 计算差异的铁死亡相关基因 
################################################################
# step2.3.1: 自己计算差异的铁死亡相关基因 
# 阙值选择；|Fold Change|>=2 & p <= 0.05
################################################################
tumor.sampleID = matadata$sampleID[matadata$group=="tumor"]
normal.sampleID = matadata$sampleID[matadata$group!="tumor"]
tumor.cpm = 2^exprMat[,tumor.sampleID]
normal.cpm = 2^exprMat[,normal.sampleID]
log2FC = log2(rowMeans(tumor.cpm)/rowMeans(normal.cpm))

# x = "ENSG00000001084"
# result = {}
# for(x in genelist){
#   temp = wilcox.test(x= tumor.cpm[x,],y=normal.cpm[x,])
#   result[[x]] = temp$p.value
# }
pvalue = unlist(lapply(genelist,function(x) {
  temp = wilcox.test(tumor.cpm[x,],normal.cpm[x,])
  return(temp$p.value)
}))

names(pvalue)  = genelist

diff_gene = data.frame("log2FC" = log2FC,"pvalue" = pvalue[names(log2FC)],row.names = names(log2FC))
diff_gene$FDR = p.adjust(diff_gene$pvalue,method = "BH") #fdr

write.table(diff_gene,"diff_ferroptosis_gene.xls",sep = "\t")

################################################################
# step2.3.1: 利用DEGSeq2计算差异的铁死亡相关基因 
diff.gene = row.names(diff_gene)[diff_gene$pvalue<=0.05]
#################################################
# tumor与normal之间差异基因与铁死亡相关基因做Venn图
genelist = diff.gene
logcpm = logcpm[,colnames(logcpm) %in% matadata$sampleID[matadata$group=="tumor"] ]
exprMat=logcpm[rownames(logcpm) %in% genelist,]
#################################################
# step2.3.2: 无监督聚类 
#################################################
library(ConsensusClusterPlus)
#为了和临床数据匹配，去掉病例后的01/06
exprMat = exprMat[,!duplicated(gsub("-01|-06","",colnames(exprMat)))]
colnames(exprMat) = gsub("-01|-06","",colnames(exprMat))
#读取临床数据
clinical = read.delim("rawdata/BRAC/Merge_clinical.txt",row.names = 1,check.names = FALSE)
#读取胰腺癌分型数据
clinical_type = read.delim("rawData/BRAC/infomation.subtype.xls",row.names = 1,check.names = FALSE)
clinical = clinical[rownames(clinical) %in% clinical_type$patient,]
clinical_type = clinical_type[clinical_type$patient %in% rownames(clinical),]
type = clinical_type$BRCA_Subtype_PAM50
names(type) = clinical_type$patient
type = data.frame(type)
clinical = merge(clinical,type,by = "row.names",all= F)
rownames(clinical)= clinical$Row.names
#将临床数据的样本于表达样本一至
clinical = clinical[rownames(clinical) %in% colnames(exprMat),]
exprMat = exprMat[,colnames(exprMat) %in% rownames(clinical)]
##################################################################
#无监督聚类
sHc<-hclust(ddist<-dist(t(exprMat)),method = "ward.D2")
res<-ConsensusClusterPlus(ddist,maxK=10,pItem=0.98,
                          reps=1000,title="me_consensus_k7_1000",clusterAlg="hc",
                          innerLinkage="ward.D2",finalLinkage="ward.D2",
                          plot="pdf",writeTable=TRUE)

#以下是构建每个每个选项
#1.1 TME3.subtype
TME3.subtype = sapply(res[[3]]$consensusClass,function(x)paste("Cluster",x,sep=""))
#1.2 生存状态
Fustat =clinical$vital_status
names(Fustat)=rownames(clinical)
#1.3年龄
clinical$A17_Age
clinical$A17_Age2= ifelse(clinical$A17_Age<=40,"20-40",clinical$A17_Age)
clinical$A17_Age2 = ifelse(clinical$A17_Age>40 &clinical$A17_Age<=60,"40-60",clinical$A17_Age2)
clinical$A17_Age2 = ifelse(clinical$A17_Age>60 &clinical$A17_Age<=80,"60-80",clinical$A17_Age2)
clinical$A17_Age2 = ifelse(clinical$A17_Age>80 &clinical$A17_Age<=100,"80-100",clinical$A17_Age2)
clinical$A17_Age2
Age = clinical$A17_Age2
names(Age)=rownames(clinical)
#1.4淋巴结转移
N = clinical$A4_N
N = gsub("N3c","N3",N)
N = gsub("N3b","N3",N)
N = gsub("N3a","N3",N)
N = gsub("N2c","N2",N)
N = gsub("N2b","N2",N)
N = gsub("N2a","N2",N)
N = gsub("N1mi","N1",N)
N = gsub("N1c","N1",N)
N = gsub("N1b","N1",N)
N = gsub("N1a","N1",N)
N = gsub("NX","N3",N)
names(N)=rownames(clinical)
#1.5 肿瘤进展情况
stage = clinical$A6_Stage
stage = gsub("Stage X","Stage IV",stage)
stage= gsub("Stage IIIA","Stage III",stage)
stage= gsub("Stage IIIB","Stage III",stage)
stage= gsub("Stage IIIC","Stage III",stage)
stage= gsub("Stage IIA","Stage II",stage)
stage= gsub("Stage IIB","Stage II",stage)
stage= gsub("Stage IIC","Stage II",stage)
stage= gsub("Stage IA","Stage I",stage)
stage= gsub("Stage IB","Stage I",stage)
stage= gsub("Stage IC","Stage I",stage)
names(stage)=rownames(clinical)
#1.6 肿瘤分型情况
type = clinical$type
names(type)=rownames(clinical)
#构建data.frame
annotation.col=data.frame(TME3.subtype,Fustat,Age,stage,N,type)
write.table(annotation.col,"annotation.col_type.xls",sep = "\t")

#######################################################################
# step2.3.3: 构建heatmap图 
#plot
htmat=scale_mat(exprMat,"row")
#每组的颜色可以不用设置，默认设置如下,也可以自己添加
annotation.col = annotation.col[order(annotation.col$TME3.subtype),]
col_ha=HeatmapAnnotation(df=annotation.col[,1:6],
                         col=list(
                           Fustat=c("Alive"="#D65F4C","Dead"="#689e46"),
                           TME3.subtype=c("Cluster1"="red","Cluster2"="blue","Cluster3"="green"),
                           stage = c("Stage I"="#EEDFCC","Stage II"="#66CDAA","Stage III"="#CD661D","Stage IV"="#8B0A50"),
                           Age   = c("20-40"="#009100","40-60"="#007979","60-80"="#004B97","80-100"="#0000C6"),
                           type  = c("Basal"="#5151A2","Her2"="#FF8040","LumA"="#613030","LumB"="#00E3E3","Normal"="#00A600"),
                           N  = c("N0"="#009100","N0"="#00E3E3","N1"="#613030","N2"="#0072E3","N3"="#5B00AE")),
                         na_col="lightgray",
                         annotation_name_gp=gpar(fontsize=6),
                         simple_anno_size=unit(3,"mm"),
                         annotation_legend_param=list(title_gp=gpar(fontsize=6,fontface="bold"),
                                                      labels_gp=gpar(fontsize=6),direction="horizontal")
)
p =Heatmap(htmat[,row.names(annotation.col)[order(annotation.col$TME3.subtype)]],
           top_annotation = col_ha,
           # row_split = annotation.row[row.names(htmat),1],
           cluster_rows = T,cluster_columns = F,
           height = unit(6, "cm"),
           # labels = c("stromal","adaptive", "innate"),
           # labels_gp = gpar(col = "white", fontsize = 10))),
           row_names_gp = gpar(fontsize = 6),show_column_names = F,
           show_parent_dend_line = FALSE,
           row_title = NULL,
           col = circlize::colorRamp2(c(-2, 0, 2), c("navyblue", "white",
                                                     "firebrick")),
           name = "NES"
)
p
###############################################################################
