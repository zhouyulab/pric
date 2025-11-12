##########################################
rm(list = ls())
library(readr)
library(dplyr)
library(DESeq2)
library(cowplot)
library(ggplot2)
library(argparse)
library(ggsignif)

args <- commandArgs(TRUE)
f_Control_Rep1 <- args[1]
f_Control_Rep2 <- args[2]
f_VIGS_Rep1 <- args[3]
f_VIGS_Rep2 <- args[4]
outtsv <- args[5]
outpdf <- args[6]

f_Control_Rep1 <- "./VIGS_RNA/results/featureCount/Control_Rep1.txt"
f_Control_Rep2 <- "./VIGS_RNA/results/featureCount/Control_Rep2.txt"
f_VIGS_Rep1 <-  "./VIGS_RNA/results/featureCount/VIGS3_Rep1.txt"
f_VIGS_Rep2 <-  "./VIGS_RNA/results/featureCount/VIGS3_Rep2.txt"
outtsv <- "./VIGS_RNA/results/featureCount/VIGS3.DEseq.txt"
outpdf <- "./VIGS_RNA/results/featureCount/VIGS3.DEseq.pdf"


df.Control.Rep1 = read_delim(f_Control_Rep1, "\t", skip = 1) 
colnames(df.Control.Rep1)[7] = "Control_Rep1"

df.Control.Rep2 = read_delim(f_Control_Rep2, "\t", skip = 1) 
colnames(df.Control.Rep2)[7] = "Control_Rep2"

df.VIGS.Rep1 = read_delim(f_VIGS_Rep1, "\t", skip = 1) 
colnames(df.VIGS.Rep1)[7] = "VIGS_Rep1"

df.VIGS.Rep2 = read_delim(f_VIGS_Rep2, "\t", skip = 1) 
colnames(df.VIGS.Rep2)[7] = "VIGS_Rep2"

df = Reduce(inner_join, list(df.Control.Rep1, df.Control.Rep2, df.VIGS.Rep1, df.VIGS.Rep2))
deseq.count = df[,c("Control_Rep1","Control_Rep2","VIGS_Rep1", "VIGS_Rep2","Geneid")]
deseq.count = deseq.count[rowSums(deseq.count[,c(1:4)])>100,]
deseq.count <- deseq.count[rowSums(deseq.count[,c(1:4)]==0) ==0,]
rowname <- deseq.count$Geneid
deseq.count <- deseq.count[,c(1:4)]
rownames(deseq.count) <- rowname
design_data <- data.frame(row.names = colnames(deseq.count), condition = c("Control", "Control","VIGS", "VIGS"))
dds <- DESeqDataSetFromMatrix(countData = deseq.count, colData = design_data, design= ~condition)
dds_rees <- DESeq(dds)

res <- results(dds_rees, contrast = c('condition', 'VIGS', 'Control'))
res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res <- res[!(is.na(res$padj)), ]
res$padj[res$padj<= 0 ]<-  min(res$padj[res$padj>0])
res$change = "nochange"
res$change[res$padj<0.05 & res$log2FoldChange> 1] <- "up"
res$change[res$padj<0.05 & res$log2FoldChange < -1] <- "down"
res$GeneName <- rownames(res)
write_tsv(res, outtsv)


df = res
df_count = df %>% group_by(change) %>% dplyr::summarise(n=n()) %>% transform(label = sprintf("%s (%s)", change, n))
p <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(pvalue), color=change))+
  geom_point(size=0.8)+
  scale_color_manual(values = c("#3952a3", "grey80", "#ef2224"), label = df_count$label)+
  labs(x= "expression foldchange (log2)", y="log10(p-value)")+
  scale_x_continuous(limits = c(-15, 15))+
  scale_x_continuous(limits = c(-15, 15))+
  geom_vline(xintercept = c(-1, 1), linetype="dashed",color="grey50")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed",color="grey50")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = c(.2, .8))
p
ggsave(outpdf, p, width=4, height =3)






###################################
rm(list = ls())
library(readr)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(GenomicFeatures)


f_Control_Rep1 <- "./VIGS_RNA/results/featureCount/Control_Rep1.txt"
f_Control_Rep2 <- "./VIGS_RNA/results/featureCount/Control_Rep2.txt"
f_VIGS1_Rep1 <-  "./VIGS_RNA/results/featureCount/VIGS1_Rep1.txt"
f_VIGS1_Rep2 <-  "./VIGS_RNA/results/featureCount/VIGS1_Rep2.txt"
f_VIGS2_Rep1 <-  "./VIGS_RNA/results/featureCount/VIGS2_Rep1.txt"
f_VIGS2_Rep2 <-  "./VIGS_RNA/results/featureCount/VIGS2_Rep2.txt"

f_gtf = "./Supplymental/Ghir_genome/Ghirsutum_HAU_gene_model.gtf"
f_outtxt = "./VIGS_RNA/results/featureCount/DESeq2.VIGS2.RNA.FPKM.txt"
f_out_normal = "./VIGS_RNA/results/featureCount/DESeq2.VIGS2.RNA.normal.count.txt"

txdb <- makeTxDbFromGFF(f_gtf, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
exon_gene_sizes <- sum(width(reduce(ebg)))

df_Control_rep1 <-  read_delim(f_Control_Rep1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_Control_rep1)[7] <- "Control_Rep1"
df_Control_rep2 <-  read_delim(f_Control_Rep2, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_Control_rep2)[7] <- "Control_Rep2"

df_VIGS1_rep1 <-  read_delim(f_VIGS1_Rep1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_VIGS1_rep1)[7] <- "VIGS1_Rep1"
df_VIGS1_rep2 <-  read_delim(f_VIGS1_Rep2, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_VIGS1_rep2)[7] <- "VIGS1_Rep2"

df_VIGS2_rep1 <-  read_delim(f_VIGS2_Rep1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_VIGS2_rep1)[7] <- "VIGS2_Rep1"
df_VIGS2_rep2 <-  read_delim(f_VIGS2_Rep2, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
colnames(df_VIGS2_rep2)[7] <- "VIGS2_Rep2"

df = Reduce(inner_join, list(df_Control_rep1, df_Control_rep2, df_VIGS1_rep1, df_VIGS1_rep2, df_VIGS2_rep1, df_VIGS2_rep2))
deseq.count = df[,c("Control_Rep1","Control_Rep2","VIGS1_Rep1", "VIGS1_Rep2","VIGS2_Rep1", "VIGS2_Rep2","Geneid")]
deseq.count = deseq.count[rowSums(deseq.count[,c(1:6)])>100,]
deseq.count <- deseq.count[rowSums(deseq.count[,c(1:6)]==0) ==0,]
rowname <- deseq.count$Geneid
deseq.count <- deseq.count[,c(1:6)]
rownames(deseq.count) <- rowname
design_data <- data.frame(row.names = colnames(deseq.count), condition = c("Control", "Control","VIGS", "VIGS", "VIGS", "VIGS"))
dds <- DESeqDataSetFromMatrix(countData = deseq.count, colData = design_data, design= ~condition)
dds_rees <- DESeq(dds)

exon_gene_sizes = exon_gene_sizes[rownames(deseq.count)]
mcols(dds)$basepairs <- exon_gene_sizes
df.fpkm = fpkm(dds)
df.fpkm = as.data.frame(df.fpkm)
df.fpkm$GeneName = rownames(df.fpkm)

