rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

load_df <- function(f_in, sepcies, samp){
  df = read_delim(f_in, "\t")
  df$zcore = (df$Reads - mean(df$Reads))/var(df$Reads)
  df = df %>% transform(Sepcies=sepcies, Sample=samp)
  return(df)
}


f_AD_WT_H3K27ac = "./pRIC/results/analysis/AD_WT_merge/AD_WT_merge.Intergenic.H3K27ac.txt"
f_AA_WT_H3K27ac = "./pRIC/results/analysis/A_WT_merge.Intergenic.H3K27ac.txt"
f_AD_WT_H3K4me3 = "./pRIC/results/analysis/AD_WT_merge/AD_WT_merge.Intergenic.H3K4me3.txt"
f_AA_WT_H3K4me3 = "./pRIC/results/analysis/A_WT_merge.Intergenic.H3K4me3.txt"

f_outtxt = "./figures/results/Figure3/Figure3A/Figure3A.txt"
f_outpdf = "./figures/results/Figure3/Figure3A/Figure3A.pdf"

df.AD.WT.H3K27ac = load_df(f_AD_WT_H3K27ac, "tetraploid", "WT") 
df.AA.WT.H3K27ac = load_df(f_AA_WT_H3K27ac, "diploid", "WT")

df.AD.WT.H3K4me3 = load_df(f_AD_WT_H3K4me3, "tetraploid", "WT")
df.AA.WT.H3K4me3 = load_df(f_AA_WT_H3K4me3, "diploid", "WT")

df = Reduce(rbind, list(df.AD.WT.H3K27ac, df.AA.WT.H3K27ac, df.AD.WT.H3K4me3, df.AA.WT.H3K4me3))
df = df %>% transform(label=paste(Sepcies, Sample, sep = "_"))

p <- ggplot(data = df, aes(x=Index, y=zcore, color=Sepcies))+
  facet_wrap(vars(Peak))+
  geom_smooth()+
  scale_x_continuous(breaks = c(0,1000,2000), labels = c("-1kb", "peak center", "1kb"))+
  labs(x="", y="enrichment of chimeric-reads \n in intergenic region (z-score)")+
  scale_color_manual(values = c("#DC3B34FF", "#04225CFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"))
p

ggsave(f_outpdf, p, width = 5, height = 2)
write_tsv(df, f_outtxt)

