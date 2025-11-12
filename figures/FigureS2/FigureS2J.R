rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggplot2)
library(cowplot)

f_hic = "/data2/zhoulab/wenmiaomiao/cotton/HiC-seq/results/analysis/HiC_gene_Chimeric.txt"
f_RIC = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge.intra.Chimeric.txt"
f_outtxt = "./figures/results/FigureS2/FigureS2J/FigureS2J.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2J/FigureS2J.pdf"


hic_total = 33976439
ric_total = 21737262

df.hic = read_delim(f_hic, "\t") 
df.hic = df.hic %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% transform(Sample="HiC-seq", total=hic_total)

df.RIC = read_delim(f_RIC, "\t") 
df.RIC = df.RIC %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% transform(Sample="RIC-seq", total = ric_total)

df = rbind(df.hic, df.RIC) %>% filter(GeneType %in% c("protein", "lncRNA", "ncRNA"))
df = df %>% transform(NormReads = ReadNum*100000000/total)
df_count = df %>% group_by(GeneType, Sample) %>% summarise(n=n())

p1 <- ggplot(data = df[df$GeneType=="lncRNA",], aes(x=Sample, y=log2(NormReads), fill=Sample))+
  facet_wrap(vars(GeneType), scales = "free")+
  geom_boxplot(aes(fill=Sample), width=0.6, outlier.shape = NA, notch=FALSE, linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill=Sample), width=0.6, outlier.shape = NA, notch=FALSE)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,color="black")+
  scale_fill_manual(values = c("grey70", "#EE2C2C"))+
  geom_text(data = df_count[df_count$GeneType=="lncRNA",], aes(x=Sample, y=-1, label=n))+
  geom_signif(comparisons = list(c("HiC-seq", "RIC-seq")), test = wilcox.test, y_position = c(12, 13, 8), 
              map_signif_level = T, color = "black",textsize = 4,
              tip_length=0.005, step_increase=0.001, margin_top = 0.05) + 
  labs(x="", y="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p1

p2 <- ggplot(data = df[df$GeneType=="protein",], aes(x=Sample, y=log2(NormReads), fill=Sample))+
  facet_wrap(vars(GeneType), scales = "free")+
  geom_boxplot(aes(fill=Sample), width=0.6, outlier.shape = NA, notch=FALSE, linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill=Sample), width=0.6, outlier.shape = NA, notch=FALSE)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,color="black")+
  scale_fill_manual(values = c("grey70", "#EE2C2C"))+
  geom_text(data = df_count[df_count$GeneType=="protein",], aes(x=Sample, y=-1, label=n))+
  geom_signif(comparisons = list(c("HiC-seq", "RIC-seq")), test = wilcox.test, y_position = c(12.5), 
              map_signif_level = T, color = "black",textsize = 4,
              tip_length=0.005, step_increase=0.001, margin_top = 0.05) + 
  labs(x="", y="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p2


p <- plot_grid(p1,p2, ncol = 2)
p
ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
