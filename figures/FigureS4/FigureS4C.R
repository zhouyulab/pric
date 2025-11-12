rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(ggsignif)

f_in = "./ChIP-seq/results/feature/enhancer.reads.norm.txt"
f_outtxt = "./figures/results/FigureS4/FigureS4C/Figure4C.txt"
f_outpdf = "./figures/results/FigureS4/FigureS4C/Figure4C.pdf"

df = read_delim(f_in, "\t")
df = df %>% transform(WT = (RNA_WT_Rep1+RNA_WT_Rep2)/2, fl = (RNA_fl_Rep1+RNA_fl_Rep2)/2, PRO_WT = (PRO_WT_Rep1+PRO_WT_Rep2)/2)
df = df %>% transform(Group = ifelse(log10(WT+1)<log10(2), "Unstable", "Stable"))
df.count = df %>% group_by(Group) %>% summarise(num=n()) %>% transform(Total = nrow(df)) %>% transform(percent=num/Total) %>%
  transform(label=paste(Group, "(",round(percent*100, 2), "%", ")", sep = ""))
df = inner_join(df, df.count)

p <- ggplot(data = df)+
  geom_histogram(aes(x=log10(WT+1), fill=label), bins = 30, color="black")+
  scale_y_continuous(expand = c(0,0))+
  geom_vline(xintercept = c(log10(2)), linetype="dashed")+
  scale_fill_manual(values = c("#CCCCFF","#FFCCCC"))+
  labs(x="RNA-seq reads per enhancer (log10 + 1)", y="Counts")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(colour = "black"))
p

ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
