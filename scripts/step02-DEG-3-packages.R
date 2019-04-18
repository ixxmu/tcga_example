library("BiocParallel")
register(MulticoreParam(2)) ##我的是两个核心 ，貌似有个检测核心个数找不找代码了…………
rm(list=ls())
options(stringsAsFactors = F)
library(DESeq2)
library(stringr)
getwd='../Rdata/'
Figure_dir='../figures/'
# 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
load( file = 
        file.path(getwd,'TCGA-KIRC-miRNA-example.Rdata')
)

expr_raw <- miRNA_tcga_xena # miRNA信息赋值给主变量，
## ** 主变量：为整个代买跑下来几乎约定俗成的变量，
rownames(expr_raw)<-expr_raw[,1]
expr_raw <- expr_raw[,-1]
meta <- miRNA_clinical # 临床信息提取
rownames(meta)<-str_replace_all(meta[,1],"-","." )## str_replace_all 替换函数
meta <-meta[,-1]
meta <- meta[str_sub(colnames(expr_raw),1,15),]##从meta里面提取所需要的样本
dim(expr)
dim(meta)
expr <- 2^expr_raw-1
exprSet<- expr##  # 为了DESeq2使用的数据；
exprSet<-round(exprSet)  ## 或者使用另一个函数，expr <- ceiling(expr)## # 取整数 为了DESeq2分析；
## gdc tcga 里面每一项都会标注数据是如何存储的，查看只有转换为原始counts进行后续计算；

# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)
group_list <- factor(group_list)
exprSet=na.omit(expr)
source('~/r_prac/pre_dara/functions1.R')

### ---------------
###
### Firstly run DESeq2 
###
### ---------------

if(T){
  library(DESeq2)
  
  (colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  tmp_f=file.path(getwd(),'TCGA-KIRC-miRNA-DESeq2-dds.Rdata')
  if(!file.exists(tmp_f)){
    dds <- DESeq(dds)
    save(dds,file = tmp_f)
  }
  load(file = tmp_f)
  res <- results(dds, 
                 contrast=c("group_list","tumor","normal"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  
  nrDEG=DESeq2_DEG[,c(2,6)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')  
  draw_h_v(exprSet,nrDEG,'DEseq2',group_list,1)
}

### ---------------
###
### Then run edgeR 
###
### ---------------
if(T){
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG =nrDEG 
  nrDEG=edgeR_DEG[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'edgeR',group_list,1)
  
}


### ---------------
###
### Lastly run voom from limma
###
### --------------- 
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('tumor-normal'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='tumor-normal', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'limma',group_list,1)
  
}
save(DEG,DEG_limma_voom,edgeR_DEG,group_list,expr,expr_raw,exprSet,meta,nrDEG,nrDEG1,file="TCGA-KIRC-miRNA-DEG_results2.Rdata") #保存标量吧
tmp_f=file.path(getwd(),'TCGA-KIRC-miRNA-DEG_results.Rdata')

if(file.exists(tmp_f)){
  save(DEG_limma_voom,DESeq2_DEG,edgeR_DEG, file = tmp_f)
  
}else{
  load(file = tmp_f) 
}



nrDEG1=DEG_limma_voom[,c(1,4)]
colnames(nrDEG1)=c('log2FoldChange','pvalue') 

nrDEG2=edgeR_DEG[,c(1,5)]
colnames(nrDEG2)=c('log2FoldChange','pvalue') 

nrDEG3=DESeq2_DEG[,c(2,6)]
colnames(nrDEG3)=c('log2FoldChange','pvalue')  

mi=unique(c(rownames(nrDEG1),rownames(nrDEG1),rownames(nrDEG1)))
lf=data.frame(lf1=nrDEG1[mi,1],
              lf2=nrDEG2[mi,1],
              lf3=nrDEG3[mi,1])
cor(na.omit(lf)) 

# 可以看到采取不同R包，会有不同的归一化算法，这样算到的logFC会稍微有差异。


