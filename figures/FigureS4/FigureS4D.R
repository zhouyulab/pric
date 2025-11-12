rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(reshape2)
library(ggsignif)

f_in = "./AD_ChIP/results/feature/enhancer.reads.norm.txt"
f_outtxt = "./figures/results/FigureS4/FigureS4D/Figure4D.txt"
f_outpdf = "./figures/results/FigureS4/FigureS4D/Figure4D.pdf"

df = read_delim(f_in, "\t")
df = df %>% transform(WT = (RNA_WT_Rep1+RNA_WT_Rep2)/2, fl = (RNA_fl_Rep1+RNA_fl_Rep2)/2, PRO_WT = (PRO_WT_Rep1+PRO_WT_Rep2)/2)
df = df %>% transform(Group = ifelse(log10(WT+1)<log10(2), "Unstable", "Stable"))
f_signal = "./AD_ChIP/results/feature/enhancer.feature.txt"
df.signal = read_delim(f_signal, "\t") %>% select(c("Enhancer", "H3K4me4", "H3K27ac", "ATAC", "DNase"))

df.merge = inner_join(df, df.signal)
df.reshape = melt(df.merge, measure.vars = c("H3K4me4", "H3K27ac", "DNase"))
df.reshape$value = as.numeric(df.reshape$value )
df.reshape = na.omit(df.reshape)

p <- ggplot(data = df.reshape, aes(x=Group, y=log2(value)))+
  facet_wrap(vars(variable), scales = "free", ncol = 4)+
  geom_boxplot(aes(fill=Group), width=0.6, outlier.shape = NA, notch=FALSE, linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill=Group), width=0.6, outlier.shape = NA, notch=FALSE)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,color="black")+
  # labs(y="Normalized PRO-CAP reads")+
  scale_fill_manual(values = c("#CCCCFF","#FFCCCC"))+
  geom_signif(comparisons = list(c("Stable", "Unstable")),  test = wilcox.test, y_position = c(0),
              map_signif_level = F, color = "black",textsize = 4,
              tip_length=0.01, step_increase=0.004, margin_top = 0.05) +
  theme_bw()+
  labs(x="")+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        legend.position = "none")
p

ggsave(f_outpdf, p, width = 4.2, height = 2.5)
write_tsv(df.reshape, f_outtxt)
