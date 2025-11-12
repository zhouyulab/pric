rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggpointdensity)

load_inter <- function(f_inter, Rep){
  df = read_delim(f_inter, "\t")
  colnames(df) <- c("GeneName1", "GeneName2", Rep)
  return(df)
}

load_intral <- function(f_intral, Rep){
  df = read_delim(f_intral, "\t")
  colnames(df) <- c("GeneName", Rep)
  return(df)
}

plot_cor <- function(df, samp, group){
  value1 = log2(df$Rep1)
  value2 = log2(df$Rep2)
  cor = cor(df$Rep1, df$Rep2)
  p <- ggplot(data = df, aes(x=log2(Rep1), y=log2(Rep2)))+
    geom_pointdensity(size=0.8)+
    ggtitle(label = sprintf("%s-%s", samp, group))+
    labs(x= sprintf("chimeric-reads of Rep1 in %s", samp), y=sprintf("chimeric-reads of Rep2 in %s", samp))+
    annotate("text", x=2, y=9, label=sprintf("R = %s", round(cor, 4)))+
    scale_color_gradient(low="orange", high="white")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          legend.position = "none")
  return(p)
}

f_rep1_inter = "./pRIC/results/analysis/A_WT_Rep1.inter.Chimeric.txt"
f_rep1_intral = "./pRIC/results/analysis/A_WT_Rep1.intra.Chimeric.txt"
f_rep2_inter = "./pRIC/results/analysis/A_WT_Rep2.inter.Chimeric.txt"
f_rep2_intral = "./pRIC/results/analysis/A_WT_Rep2.intra.Chimeric.txt"

f_outtxt = "./figures/results/FigureS2/FigureS2B/FigureS2B.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2B/FigureS2B.pdf"


df.rep1.inter = load_inter(f_rep1_inter, "Rep1")
df.rep2.inter = load_inter(f_rep2_inter, "Rep2")
df.inter = inner_join(df.rep1.inter, df.rep2.inter)
p_inter = plot_cor(df = df.inter, "AA", "inter")
p_inter

df.rep1.intral = load_intral(f_rep1_intral, "Rep1")
df.rep2.intral = load_intral(f_rep2_intral, "Rep2")
df.intral = inner_join(df.rep1.intral, df.rep2.intral)
p_intral = plot_cor(df = df.intral, "AA", "intra")
p_intral

df = df.intral %>% transform(GeneName2=GeneName) %>% dplyr::rename("GeneName1"="GeneName") %>% rbind(df.inter)
p = plot_cor(df = df, "AA", "all")
p

p_merge = plot_grid(p_inter, p_intral, p, ncol = 3)
p_merge

ggsave(f_outpdf, p_merge, width = 6, height = 2.2)
write_tsv(df, f_outtxt)
