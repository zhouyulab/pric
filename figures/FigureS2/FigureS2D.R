rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggblanket)
library(ggdensity)
library(ggpointdensity)

plot_cor <- function(df){
  value1 = log2(df$ReadNum)
  value2 = log2(df$expression)
  cor_value = round(cor(value1, value2), 4)
  p <- ggplot(data = df, aes(x=log2(ReadNum), y=log2(expression)))+
    geom_pointdensity(size=0.8)+
    geom_hdr_lines(probs = c(0.99, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1), linewidth=0.4, color="black", n=50) +
    annotate("text", x=1, y=9, label=sprintf("R = %s", cor_value), hjust=0)+
    scale_color_gradient(low = "#FFA500", high ="white")+
    labs(x="Chimeric reads (log2)", y="Gene expression (log2)")+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"))
  p
  return(p)
}

f_RNA = "./AD_RNA/results/featureCount/RNA.FPKM.txt"
f_wt_chimeric = "./pRIC/results/analysis/AD_WT_merge/AD_WT_merge.intra.Chimeric.txt"
f_fl_chimeric = "./pRIC/results/analysis/AD_fl_merge/AD_fl_merge.intra.Chimeric.txt"
f_outtxt = "./figures/results/FigureS2/FigureS2D/FigureS2D.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2D/FigureS2D.pdf"

df_RNA = read_delim(f_RNA, "\t")
df_RNA = df_RNA %>% transform(WT = (WT_Rep1+WT_Rep2)/2, fl = (fl_Rep1+fl_Rep2)/2)

df.wt.chimeric = read_delim(f_wt_chimeric, "\t")
df.wt.chimeric = df.wt.chimeric %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% filter(GeneType=="protein")
df.wt.chimeric = df.wt.chimeric %>% inner_join(df_RNA[, c("GeneName", "WT")]) %>% dplyr::rename("expression"="WT")

df.fl.chimeric = read_delim(f_fl_chimeric, "\t")
df.fl.chimeric = df.fl.chimeric %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% filter(GeneType=="protein")
df.fl.chimeric = df.fl.chimeric %>% inner_join(df_RNA[, c("GeneName", "fl")]) %>% dplyr::rename("expression"="fl")

df = rbind(df.wt.chimeric, df.fl.chimeric)
p_merge = plot_cor(df)
p_merge

ggsave(f_outpdf, p_merge, width = 3.3, height = 3)
write_tsv(df, f_outtxt)


