library(ggplot2)
library(pheatmap)

treat="CTEAD"
control="CTE"

sub_gl <- subset(gl,gl$group_list %in% c(treat,control))

sub_group_list <- sub_gl$group_list
colnames(sub_gl)[1] <- 'group_list'

sub_count <-counts[,colnames(counts) %in% rownames(sub_gl)]

##筛选条件设置
log2FC_cutoff = 1
padj_cutoff = 0.05
outdir ='./'
MyDEseq(sub_count, sub_gl, sub_group_list,treat,control,log2FC_cutoff, padj_cutoff, outdir)

MyDEseq <- function(counts, gl, group_list, treat, control, log2FC_cutoff, padj_cutoff, outdir){
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = gl,
                                design = ~ group_list)
  
  dds$group_list <- relevel(dds$group_list, ref = control)   #指定 control group
  keep <- rowSums(counts(dds)) >= 1.5*ncol(counts)  #Pre-filtering ，过滤低表达基因
  dds <- dds[keep,] 
  dds <- DESeq(dds,quiet = F) 
  res <- results(dds,contrast=c("group_list", treat, control))  #指定提取为exp/ctr结果
  resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
  tempDEG <- as.data.frame(resOrdered)
  DEG_DEseq2 <- na.omit(tempDEG)
  
  #### 保存DEG数据
  save(DEG_DEseq2, file =paste0(outdir, '/',control,'_vs_',treat, '_DEG_res.Rdata')) 
  
  ## Volcano plot

  ##选取差异分析结果
  need_DEG <- DEG_DEseq2[,c(2,6)] #选取log2FoldChange, padj信息
  colnames(need_DEG) <- c('log2FoldChange','padj') 
  
  need_DEG$significance  <- as.factor(ifelse(need_DEG$padj < padj_cutoff & abs(need_DEG$log2FoldChange) > log2FC_cutoff,
                                             ifelse(need_DEG$padj < padj_cutoff & need_DEG$log2FoldChange > log2FC_cutoff ,'UP','DOWN'),'NOT'))
  
  title <- paste0(' Up :  ',nrow(need_DEG[need_DEG$significance =='UP',]) ,
                  '\n Down : ',nrow(need_DEG[need_DEG$significance =='DOWN',]),
                  '\n FoldChange >= ',round(2^log2FC_cutoff,3))
  
  pvol <- ggplot(data=need_DEG, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
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
  ggsave(pvol,filename = paste0(outdir, '/',control,'_vs_',treat, '_vol_',log2FC_cutoff,'_',padj_cutoff,'.png'),width =8,height =7.5)
  
  ## heatmap
  ##选择要展示基因表达量的数据
  
  gene_up <- rownames(need_DEG[with(need_DEG,log2FoldChange>log2FC_cutoff & padj<padj_cutoff),])
  gene_down <- rownames(need_DEG[with(need_DEG,log2FoldChange< -log2FC_cutoff & padj<padj_cutoff),])
  annotation_col = gl
  colnames(annotation_col)[1] ='group'
  
  pheatmap::pheatmap(counts[c(gene_up,gene_down),], scale="row",angle_col = 45,
                     color = colorRampPalette(c('#115699','white','#BB1E38'))(100),
                     show_colnames =T,
                     show_rownames = F,
                     fontsize = 7 ,
                     cluster_cols = F,
                     annotation_col=annotation_col,
                     filename = paste0(outdir, control,'_vs_',treat, '_heatmap_',log2FC_cutoff,'_',padj_cutoff,'.png'))

}





