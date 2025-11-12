rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(circlize)
library(cowplot)
library(argparse)

f_pCp = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000000.bin.txt"
f_id = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000000.ID.txt"
f_outtxt = "./figures/results/Figure1/Figure1D/Figure1D.txt"
f_outpdf = "./figures/results/Figure1/Figure1D/Figure1D.pdf"

df.id = read_delim(f_id, "\t")
df.id = df.id %>% group_by(Chrom) %>% summarise(Pos=max(ID))

df.pCp = read_delim(f_pCp, "\t") 
df.pCp$ReadNum[df.pCp$ReadNum>1000] = 1000
df.pCp = df.pCp %>% filter(ReadNum>1)
p <- ggplot()+
  geom_rect(data = df.pCp, aes(xmin=ID1-0.5, xmax=ID1+0.5, ymin=ID2-0.5, ymax=ID2+0.5, fill=log2(ReadNum)))+
  geom_rect(data = df.pCp,  aes(xmin=ID2-0.5, xmax=ID2+0.5, ymin=ID1-0.5, ymax=ID1+0.5, fill=log2(ReadNum)))+
  scale_fill_gradientn(colours = c("grey90", "#FFA500", "red","darkred"))+
  scale_x_continuous(expand = c(0,0), breaks = df.id$Pos, labels = df.id$Chrom, position = "top")+
  scale_y_continuous(expand = c(0,0), breaks = df.id$Pos, labels = df.id$Chrom)+
  theme_bw()+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

p
ggsave(f_outpdf, p, width = 5, height = 3.5)
write_tsv(df.pCp, f_outtxt)
