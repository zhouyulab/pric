rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

annotate_pvalue <- function(pvalue){
  if(pvalue < 0.001){
    annotation <- "***"
  }else{
    if(pvalue < 0.01){
      annotation <- "**"
    }else{
      if(pvalue < 0.05){
        annotation <- "*"
      }else{
        annotation <- "."
      }
    }
  }
  return(annotation)
}


f_reads = "./AD_RNA/results/featureCount/RNA.FPKM.txt"
f_in = "./pRIC/results/figures//EP.foldchange.tsv"
f_outtxt = "./figures/results/FigureS9/FigureS9D/FigureS9D.txt"
f_outpdf = "./figures/results/FigureS9/FigureS9D/FigureS9D.pdf"

f_gene = "./Supplymental/Ghir_genome/annotation/gene.anno.fiber.tsv"
df.gene = read_delim(f_gene, "\t")
fiber_gene = df.gene$cottonGene[df.gene$idFiber | df.gene$isC3]

df.reads = read_delim(f_reads, "\t")
df.reads = df.reads %>% transform(WT = (WT_Rep1+WT_Rep2)/2,  FL = (fl_Rep1+fl_Rep2)/2 )

df = read_delim(f_in, "\t")
df.WT = df %>% filter(Group=="WT-specific") %>% arrange(desc(log2FoldChange)) 
df.fl = df %>% filter(Group=="fl-specific") %>% arrange(log2FoldChange) 

wt_gene = df.WT$GeneName[1:30]
fl_gene = df.fl$GeneName[1:30]

df.reads.WT = df.reads %>% filter(GeneName %in% wt_gene) %>% transform(fc=WT-FL) %>% arrange(desc(fc)) %>% transform(group="WT")
df.reads.fl = df.reads %>% filter(GeneName %in% fl_gene) %>% transform(fc=WT-FL) %>% arrange(desc(fc)) %>% transform(group="fl")
df.reads.filter = rbind(df.reads.WT, df.reads.fl)
df.reads.filter$pvalue = apply(df.reads.filter[,c("WT_Rep1", "WT_Rep2", "fl_Rep1", "fl_Rep2")], 1, function(x){
  wt = as.numeric(c(x[1], x[2]))
  fl = as.numeric(c(x[3], x[4]))
  pvalue = t.test(wt,fl, var.equal = TRUE)$p.value
  pvalue_anno = annotate_pvalue(pvalue)
  return(pvalue_anno)
})

df.gene.filter = df.gene %>% filter(cottonGene %in% df.reads.filter$GeneName) %>% dplyr::select(c("cottonGene", "GeneName"))
colnames(df.gene.filter) <- c("GeneName", "Symbol")
df.reads.filter = left_join(df.reads.filter, df.gene.filter)
df.reads.filter$Symbol[is.na(df.reads.filter$Symbol)] <- "None"
df.reads.filter = df.reads.filter %>% transform(label=paste("(",Symbol, ")",GeneName, sep = ""))

dat = t(df.reads.filter[,c(1:4)])
dat = as.matrix(log2(dat))
colnames(dat) = df.reads.filter$label

col = c("higher" = "#67A9CFFF", "lower" = "#EF8A62FF")
top_anno = columnAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                             labels = c("higher", "lower"), 
                                             panel_fun = function(index, nm){grid.rect(gp = gpar(fill = col[nm]))
                                               grid.text(nm, 0.5, 0.5)},
                                             labels_gp = gpar(col = "white", fontsize = 10, width=1),
                                             height = unit(0.5, "cm")),
                            pvalue = anno_simple(df.reads.filter$pvalue,
                                                 pch = df.reads.filter$pvalue,
                                                 border = FALSE,
                                                 na_col = "grey")                
)


pdf(f_outpdf, width = 12, height = 5)
Heatmap(dat,
        name = "log2(expression)",
        border_gp = gpar(col = "black"),
        rect_gp = gpar(col = "black"),
        column_split = df.reads.filter$group,
        row_names_side = "left",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_names_centered = TRUE,
        cluster_rows = TRUE,
        col = colorRamp2(c(-5, 2, 10), c("#72A7CF", "white", "#FF7256")),
        top_annotation = top_anno,
        width = unit(24, "cm"),
        height = unit(5, "cm"))


dev.off()
sum(df.reads.WT$GeneName %in% fiber_gene)

outdf = df.reads.filter %>% dplyr::select(c("GeneName", "WT", "FL", "group", "pvalue", "Symbol"))
write_tsv(outdf, f_outtxt)
