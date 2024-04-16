setwd('D:\\05_project\\dem_tbi\\CTE_AD\\Dataset_1')

#### preprocess ####
###环境设置
rm(list=ls())
options(stringsAsFactors = F) 
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件

## 对counts进行处理筛选得到表达矩阵 
a1 <- fread('./all_2.id.txt',
            header = T,data.table = F)#载入counts，第一列设置为列名
colnames(a1)
counts <- a1[,7:ncol(a1)] #截取样本基因表达量的counts部分作为counts
rownames(counts) <- a1$Geneid #将基因名作为行名
#更改样品名
colnames(counts)
colnames(counts) <- gsub('/home/test/align/bam/','', #删除样品名前缀
                         gsub('.bam','',  colnames(counts))) #删除样品名后缀
colnames(counts)


## 导入或构建样本信息,  进行列样品名的重命名和分组
# b <- read.csv('./SraAccList.csv')
# b
# name_list <- b$source_name[match(colnames(counts),b$Run)]; name_list  #选择所需要的样品信息列
# nlgl <- data.frame(row.names=colnames(counts),
#                    name_list=name_list,
#                    group_list=gsub("_.*", "", name_list))
# 
# name_list <- nlgl$name_list
# colnames(counts) <- name_list #更改样品名
# group_list <- nlgl$group_list

gl <- data.frame(row.names=colnames(counts), #构建样品名与分组对应的数据框
                 group_list=c(rep("AD",10),rep("CTE",8),rep("CTEAD",6),rep("NORMAL",10)))

group_list <- gl$group_list


## counts，TPM转化 
# # 注意需要转化的是未经筛选的counts原始矩阵
# ### 从featurecounts 原始输出文件counts.txt中提取Geneid、Length(转录本长度)，计算tpm
# geneid_efflen <- subset(a1,select = c("Geneid","Length"))
# colnames(geneid_efflen) <- c("geneid","efflen")  
# 
# ### 取出counts中geneid对应的efflen
# efflen <- geneid_efflen[match(rownames(counts),
#                               geneid_efflen$geneid),
#                         "efflen"]
# 
# ### 计算 TPM 公式
# #TPM (Transcripts Per Kilobase Million)  每千个碱基的转录每百万映射读取的Transcripts
# counts2TPM <- function(count=count, efflength=efflen){
#   RPK <- count/(efflength/1000)   #每千碱基reads (Reads Per Kilobase) 长度标准化
#   PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
#   RPK/PMSC_rpk              
# }  
# 
# tpm <- as.data.frame(apply(counts,2,counts2TPM))
# colSums(tpm)

#合并所有重复symbol
g2s <- fread('./g2s_Homo_sapiens.GRCh38.111.txt',header = F,data.table = F) #载入从gencode的gtf文件中提取的信息文件
colnames(g2s) <- c("geneid","symbol")

symbol <- g2s[match(rownames(counts),g2s$geneid),"symbol"] #匹配counts行名对应的symbol
table(duplicated(symbol))  #统计重复基因名

# #### id ####
# library(biomaRt)
# 
# ensembl_id <- rownames(counts)
# mart <- useMart("ensembl","hsapiens_gene_ensembl")
# dataset = listDatasets(mart)
# mydataset = useDataset("hsapiens_gene_ensembl",mart = mart)
# symbols <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
#                     filters = 'ensembl_gene_id', values = ensembl_id, mart = mydataset)
# symbols <- symbols[which(symbols$external_gene_name != ""),]
# 
# countf2 <- data.frame(ensembl_gene_id = rownames(counts),counts)
# 
# countf2 <- left_join(symbols,countf2,by=c("ensembl_gene_id"="ensembl_gene_id"))
# 
# countf2 <- aggregate(x = countf2[,3:ncol(countf2)],   #此时exprSet的第三列开始是表达矩阵内容
#                      by = list(symbol = countf2$external_gene_name),   #按照相同symbol分组，在组内计算
#                      FUN = sum) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
#   column_to_rownames(var = 'symbol')
# write.csv(countf2, './count_mat_1.csv')
#### id ####


###使用aggregate根据symbol列中的相同基因进行合并 

counts <- aggregate(counts, by=list(symbol), FUN=sum)

counts <- counts[counts$Group.1 != "", ]

rownames(counts) <- counts[,1]

counts <- counts[,-1]

# tpm <- aggregate(tpm, by=list(symbol), FUN=sum) ###使用aggregat 将symbol列中的相同基因进行合并 
# tpm <- column_to_rownames(tpm,'Group.1')

#### 初步过滤低表达基因 ####（筛选标准不唯一、依情况而定）
#筛选出至少在重复样本数量内的表达量counts大于1的行（基因）
keep_feature <- rowSums(counts>1) >= 2
table(keep_feature)  #查看筛选情况，FALSE为低表达基因数（行数），TURE为要保留基因数


# counts_filt <- counts[keep_feature, ] #替换counts为筛选后的基因矩阵（保留较高表达量的基因）
# tpm_filt <- tpm[keep_feature, ]

## 保存数据
counts_raw=counts #这里重新命名方便后续分析调用
# counts=counts_filt
# tpm=tpm_filt

save(counts_raw, counts, gl, group_list, file='./1.counts.Rdata')

#### 02_check ####
rm(list = ls())  
options(stringsAsFactors = F)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#### 载入数据 设置目录
setwd('D:\\05_project\\dem_tbi\\CTE_AD\\Dataset_1')
load(file = 'D:\\05_project\\dem_tbi\\CTE_AD\\Dataset_1/1.counts.Rdata')
dir.create("2.check")
setwd("2.check")

#### 数据预处理  # (任选以下一种作为dat即可，主要是进行样本间归一化，使得样本具有可比性)
#dat <- as.data.frame(log2(edgeR::cpm(counts)+1))  #简单归一化 CPM:Counts per million
#dat <- log2(tpm+1)
#DESeq2_normalize   rld 
group_list <- gl$group_list
if (T) { 
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = gl,
                                design = ~ group_list)
  rld <- rlog(dds, blind=FALSE)
  write.table(assay(rld),  file="Deseq2_rld.txt", sep="\t", quote=F, col.names=NA)
  dat <- as.data.frame(assay(rld))     
}

### boxplot 查看样本的基因整体表达情况
boxplot(dat,col=colors, ylab="dat", main=" normalized data ",
        outline = F, notch = F)
dev.off()

###################### hclust and Heatmap of the sample-to-sample distances ###########################
sampleDists <- dist(t(dat))   #dist默认计算矩阵行与行的距离， 因此需要转置
sampleDistMatrix <- as.matrix(sampleDists)  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  #选取热图的颜色
p0 <- pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
ggsave(p0,filename = 'check_dist.pdf',width = 7.5,height =6)
dev.off()

pdf("check_hclust.pdf")
plot(hclust(sampleDists))
dev.off()

################################# PCA检测 #####################################
#PCA查看实验和对照组情况
pca <- plotPCA(rld, ntop = 500, intgroup=c("group_list"))
ggsave(pca, filename = 'check_PCA.pdf',width = 7.5,height =6)

####################### heatmap检测——取500差异大的基因 ##########################################
cg <- names(tail(sort(apply(dat,1,sd)),500)) #取每一行的方差，从小到大排序，取最大的500个
n <- dat[cg,]
p1 <- pheatmap::pheatmap(n,show_colnames =T,show_rownames = F,
                         fontsize=7,
                         legend_breaks = -3:3,
                         #scale = "row",
                         angle_col=45,
                         annotation_col=gl) 

ggsave(p1,filename = 'check_heatmap_top500_sd.pdf',width = 7.5,height =6)
dev.off()

#######################样本相关性检测————取500高表达基因##################################
dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),]#取高表达量前500基因
M <- cor(dat_500)

p2 <-pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        fontsize=7,
                        annotation_col = gl ) 
ggsave(p2,filename = 'check_cor_top500.pdf',width = 7.5,height =6)
dev.off()

#### 03_DEG ####
rm(list = ls())  
options(stringsAsFactors = F)

setwd('D:\\05_project\\dem_tbi\\CTE_AD\\Dataset_1')
load(file = '1.counts.Rdata')
dir.create("3.DEG")
setwd("3.DEG")
library(DESeq2)
library("BiocParallel") #启用多核计算

##设定 实验组exp / 对照组ctr



exp="CTE"
ctr="NORMAL"

sub_gl <- subset(gl,gl$group_list %in% c(exp,ctr))

sub_group_list <- sub_gl$group_list
colnames(sub_gl)[1] <- 'sub_group_list'

sub_count <-counts[,colnames(counts) %in% rownames(sub_gl)]
# sub_tpm <- tpm[,colnames(tpm) %in% rownames(sub_gl)]
#### DEG using DESeq2 ####

##构建dds DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = sub_count,
                                colData = sub_gl,
                                design = ~ sub_group_list)

dds$sub_group_list <- relevel(dds$sub_group_list, ref = ctr)   #指定 control group
keep <- rowSums(counts(dds)) >= 1.5*ncol(sub_count)  #Pre-filtering ，过滤低表达基因
dds <- dds[keep,] 
dds <- DESeq(dds,quiet = F) 
res <- results(dds,contrast=c("sub_group_list", exp, ctr))  #指定提取为exp/ctr结果
resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)

#### 保存DEG数据
save(DEG_DEseq2, file ='test_DEG_results.Rdata') 

#### Volcano plot ####
library(ggplot2)
library(pheatmap)
##筛选条件设置
log2FC_cutoff = 1
padj_cutoff = 0.05
##选取差异分析结果
need_DEG <- DEG_DEseq2[,c(2,6)] #选取log2FoldChange, padj信息
colnames(need_DEG) <- c('log2FoldChange','padj') 

need_DEG$significance  <- as.factor(ifelse(need_DEG$padj < padj_cutoff & abs(need_DEG$log2FoldChange) > log2FC_cutoff,
                                           ifelse(need_DEG$padj < padj_cutoff & need_DEG$log2FoldChange > log2FC_cutoff ,'UP','DOWN'),'NOT'))

title <- paste0(' Up :  ',nrow(need_DEG[need_DEG$significance =='UP',]) ,
                '\n Down : ',nrow(need_DEG[need_DEG$significance =='DOWN',]),
                '\n FoldChange >= ',round(2^log2FC_cutoff,3))

g <- ggplot(data=need_DEG, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
  #点和背景
  geom_point(alpha=0.8, size=2) +
  theme_classic()+ #点的透明度，大小
  labs(x = 'Log2 fold change',
       y = '-Log10 Padj')+
  #标题文本
  ggtitle( title ) +
  #分区颜色                  
  scale_colour_manual(values = c('#115699','grey','#BB1E38'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(), #不显示图例标题
        legend.position="right")  #图例位置
dev.off()
ggsave(g,filename = 'volcano_padj.png',width =8,height =7.5)

#### heatmap ####
##选择要展示基因表达量的数据

gene_up <- rownames(need_DEG[with(need_DEG,log2FoldChange>log2FC_cutoff & padj<padj_cutoff),])
gene_down <- rownames(need_DEG[with(need_DEG,log2FoldChange< -log2FC_cutoff & padj<padj_cutoff),])
cg <- c(head(gene_up, 100),   #取前50 padj上下调基因名
        head(gene_down, 100))
cg <- na.omit(match(cg,rownames(dat))) 

annotation_col = sub_gl
colnames(annotation_col)[1] ='group'

pheatmap::pheatmap(sub_count[cg,], scale="row", angle_col = 45,
                   color = colorRampPalette(c('#115699','white','#BB1E38'))(100),
                   show_colnames =T,
                   show_rownames = F,
                   fontsize = 7 ,
                   cluster_cols = F,
                   annotation_col=annotation_col)


pheatmap::pheatmap(sub_count[c(gene_up,gene_down),], scale="row",angle_col = 45,
                   color = colorRampPalette(c('#115699','white','#BB1E38'))(100),
                   show_colnames =T,
                   show_rownames = F,
                   fontsize = 7 ,
                   cluster_cols = T,
                   annotation_col=annotation_col,
                   filename = 'heatmap_all_up&down_DEG_2.png')

