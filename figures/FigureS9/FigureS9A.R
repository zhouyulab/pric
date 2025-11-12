rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)
library(ggbreak)

group_chimeric <- function(readnum){
  if(readnum<=3){return("<=3")}
  if(readnum>3 & readnum<=5){return("4~5")}
  if(readnum>5 & readnum<=10){return("6~10")}
  if(readnum>10){return(">10")}
  
}

f_wt = "./pRIC/results/EP/AD_WT_merge/AD_WT_merge.EP.fisher.txt"
f_fl = "./pRIC/results/EP/AD_fl_merge/AD_fl_merge.EP.fisher.txt"
f_outtxt = "./figures/results/FigureS9/FigureS9A/FigureS9A.txt"
f_outpdf = "./figures/results/FigureS9/FigureS9A/FigureS9A.pdf"

df.wt = read_delim(f_wt, "\t") %>% transform(Sample = "WT") %>% filter(Group=="E-P")
df.wt$ChimericGroup = apply(as.data.frame(df.wt[,"Reads"]), 1, group_chimeric)
df.fl = read_delim(f_fl, "\t") %>% transform(Sample = "fl")  %>% filter(Group=="E-P")
df.fl$ChimericGroup = apply(as.data.frame(df.fl[,"Reads"]), 1, group_chimeric)

df = rbind(df.wt, df.fl)
df$ChimericGroup = factor(df$ChimericGroup, levels = c("<=3", "4~5", "6~10", ">10"))
df$Sample = factor(df$Sample, levels = c("WT", "fl"))
df = df %>% filter(ChimericGroup != "<=3")
df_count = df %>% group_by(ChimericGroup, Sample) %>% summarise(num=n())
p <- ggplot(data = df_count, aes(x=ChimericGroup, y=num, fill=Sample))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_text(aes(x=ChimericGroup, y=num+100, label=num))+
  labs(x="chimeric reads", y="E-P count") +
  # scale_y_break(breaks = c(1000, 4000))+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90),
        axis.text.y = element_text(colour = "black"),
        legend.title = element_blank())
p
ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df_count, f_outtxt)
