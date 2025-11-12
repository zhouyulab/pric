rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggbreak)

f_outtxt = "./figures/results/FigureS6/FigureS6D/FigureS6D.txt"
f_outpdf = "./figures/results/FigureS6/FigureS6D/FigureS6D.pdf"
f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.fisher.txt"
f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.EP.fisher.txt"

df.AD = read_delim(f_AD, "\t") %>% transform(Sample="Tetraploid")
df.AA = read_delim(f_AA, "\t") %>% transform(Sample="Diploid")
df = rbind(df.AD, df.AA)

df.enhancer = df %>% group_by(Sample, Enhancer) %>% summarise(num=n())
df.total =  df.enhancer  %>% group_by(Sample) %>% summarise(Total=n())
df.enhancer = df.enhancer %>% transform(Group=ifelse(num>6, ">6", num))
df.enhancer.count = df.enhancer %>% group_by(Group,Sample) %>% summarise(n=n()) %>% inner_join(df.total) %>% transform(percent=n*100/Total)

df.enhancer.count$Group = factor(df.enhancer.count$Group, levels = c("1", "2", "3", "4", "5", "6", ">6"))
df.enhancer.count$Sample = factor(df.enhancer.count$Sample, levels = c("Diploid", "Tetraploid"))
p_enhancer <- ggplot(data = df.enhancer.count, aes(x=Group, y=percent, fill=Sample))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = c("#E41A1C","#377EB8"))+
  labs(x="Number of linked promoter", y="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(.8, .6))
p_enhancer

ggsave(f_outpdf, p_enhancer, width = 3, height = 2.5)
write_tsv(df.enhancer.count, f_outtxt)
