rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggbreak)

f_outtxt = "./figures/results/FigureS6/FigureS6E/FigureS6E.txt"
f_outpdf = "./figures/results/FigureS6/FigureS6E/FigureS6E.pdf"
f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.fisher.txt"
f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.EP.fisher.txt"

df.AD = read_delim(f_AD, "\t") %>% transform(Sample="Tetraploid")
df.AA = read_delim(f_AA, "\t") %>% transform(Sample="Diploid")
df = rbind(df.AD, df.AA)

df.promoter = df %>% group_by(Sample, Promoter) %>% summarise(num=n())
df.total =  df.promoter  %>% group_by(Sample) %>% summarise(Total=n())
df.promoter = df.promoter %>% transform(Group=ifelse(num>6, ">6", num))
df.promoter.count = df.promoter %>% group_by(Group,Sample) %>% summarise(n=n()) %>% inner_join(df.total) %>% transform(percent=n*100/Total)

df.promoter.count$Group = factor(df.promoter.count$Group, levels = c("1", "2", "3", "4", "5", "6", ">6"))
df.promoter.count$Sample = factor(df.promoter.count$Sample, levels = c("Diploid", "Tetraploid"))
p_promoter <- ggplot(data = df.promoter.count, aes(x=Group, y=percent, fill=Sample))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = c("#E41A1C","#377EB8"))+
  labs(x="Number of linked enhancer", y="")+
  # scale_y_break(c(20, 60))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(.8, .6))
p_promoter

ggsave(f_outpdf, p_promoter, width = 6, height = 2.5)
write_tsv(df.promoter.count, f_outtxt)
