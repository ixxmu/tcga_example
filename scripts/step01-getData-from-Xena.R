## 
### ---------------
### ---------------

# TCGA-XENA数据库 打包：https://share.weiyun.com/56URQ3a
# TCGA-GDC-somatic  链接：https://share.weiyun.com/5fx40jk 密码：7yenp9 

rm(list=ls())
options(stringsAsFactors = F)
d='../Rdata/KIRC-UCSC-XENA/'
## too many NA in the miRNA expression matrix from XENA

if(file.exists(file.path(d,'miRNA_GA_gene.gz'))){
  
  miRNA_GA=read.table('TCGA.KIRC.sampleMap__miRNA_GA_gene.gz',header = T,sep = '\t')
  dim(miRNA_GA)
  miRNA_GA[1:4,1:4]
  na.omit(miRNA_GA)[1:4,1:4]
  dim(na.omit(miRNA_GA))
  
  miRNA_HiSeq=read.table("TCGA.KIRC.sampleMap__miRNA_HiSeq_gene.gz" ,header = T,sep = '\t')
  dim(miRNA_HiSeq)
  miRNA_HiSeq[1:4,1:4]
  na.omit(miRNA_HiSeq)[1:4,1:4]
  dim(na.omit(miRNA_HiSeq))
}

miRNA_tcga_xena=read.table("TCGA-KIRC__Xena_Matrices__TCGA-KIRC.mirna.tsv.gz" ,header = T,sep = '\t')
miRNA_clinical=read.table("TCGA.KIRC.sampleMap__KIRC_clinicalMatrix.gz" ,header = T,sep = '\t')
save(miRNA_GA,miRNA_HiSeq,miRNA_tcga_xena,miRNA_clinical,file="TCGA-KIRC-miRNA-example.Rdata")


