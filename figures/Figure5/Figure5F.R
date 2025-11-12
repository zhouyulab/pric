rm(list = ls())
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(argparse)
library(ggsignif)

load_EP_gene <- function(f_in, cutoff=0){
  EP_gene = c()
  df = read_delim(f_in, "\t") %>% filter(Group=="E-P") %>% filter(Reads>cutoff)
  for(i in 1:nrow(df)){
    gene = gsub("pro", "", df$Promoter[i])
    EP_gene = c(EP_gene, gene)
  }
  return(EP_gene)
}

f_outtxt = "./figures/results/Figure5/Figure5F/Figure5F.txt"
f_outpdf = "./figures/results/Figure5/Figure5F/Figure5F.pdf"

f_WT = "./pRIC/results/EP/AD_WT_merge/AD_WT_merge.EP.fisher.txt"
f_fl = "./pRIC/results/EP/AD_fl_merge/AD_fl_merge.EP.fisher.txt"
f_RNA = "./AD_RNA/results/featureCount/RNA.FPKM.txt"

df_RNA = read_delim(f_RNA, "\t")
df_RNA = df_RNA %>% transform(WT = (WT_Rep1+WT_Rep2)/2, fl = (fl_Rep1+fl_Rep2)/2)

WT_gene = load_EP_gene(f_WT, cutoff = 3)
fl_gene = load_EP_gene(f_fl, cutoff = 3)
common_gene = intersect(WT_gene, fl_gene)

f_DEseq = "./AD_RNA/results/featureCount/RNA.DEseq.res"
df.DEseq = read_delim(f_DEseq, "\t")
df.DEseq$Group = "None"
df.DEseq$Group[df.DEseq$GeneName %in% common_gene] <- "common"
df.DEseq$Group[df.DEseq$Group=="None" & df.DEseq$GeneName %in% WT_gene] <- "WT-specific"
df.DEseq$Group[df.DEseq$Group=="None" & df.DEseq$GeneName %in% fl_gene] <- "fl-specific"
df.DEseq = df.DEseq %>% filter(Group != "None")
df.DEseq$Group = factor(df.DEseq$Group, levels = c("WT-specific", "common", "fl-specific"))

df_count = df.DEseq %>% group_by(Group) %>% summarise(n=n())
df.point = df.DEseq %>% group_by(Group) %>% summarise(meanvalue=mean(log2FoldChange), medianvalue=median(log2FoldChange))
p <- ggplot(data = df.DEseq, aes(x=Group, y=log2FoldChange, fill=Group))+
  geom_violin()+
  geom_boxplot(width=0.3, notch=TRUE, outlier.shape = NA)+
  geom_point(data = df.point, aes(x=Group, y=medianvalue))+
  scale_fill_manual(values = c("#e8484a", "#5973a2","#5e99c7"))+
  scale_color_manual(values = c("#e8484a", "#5973a2","#5e99c7"))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x="", y="gene foldchange (WT/fl, log2)")+
  geom_hline(yintercept = 0, linetype="dashed", color="grey60")+
  geom_text(data = df_count, aes(x=Group, y=-2.5, label=n))+
  geom_text(data = df.point, aes(x=Group, y=-0.5, label=round(medianvalue, 4)))+
  geom_signif(comparisons = list(c("WT-specific", "common"), c("WT-specific", "fl-specific"), c("fl-specific", "common")),
              test = t.test, y_position = c(0.8, 1, 0.4),
              map_signif_level = F, color = "black",textsize = 4,
              tip_length=0.005, step_increase=0.015, margin_top = 0.05) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black"))
p
ggsave(f_outpdf, p, width = 2, height = 3)
write_tsv(df.DEseq, f_outtxt)
