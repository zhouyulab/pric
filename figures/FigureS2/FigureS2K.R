rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggplot2)
library(circlize)
library(cowplot)
library(argparse)
library(ComplexHeatmap)
f_pCp = "./pRIC/results/analysis/A_WT_merge/A_WT_merge_1000000.bin.txt"
f_id = "./pRIC/results/analysis/A_WT_merge/A_WT_merge_1000000.ID.txt"
f_outtxt = "./figures/results/FigureS2/FigureS2K/FigureS2K.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2K/FigureS2K.pdf"

df.id = read_delim(f_id, "\t")
df.id = df.id %>% group_by(Chrom) %>% summarise(Pos=max(ID))
df.id = df.id %>% filter(Chrom != "Chr14")

df.pCp = read_delim(f_pCp, "\t") 
df.pCp = df.pCp %>% filter(ReadNum >1)
df.pCp$ReadNum[df.pCp$ReadNum>1000] = 1000
df.pCp = df.pCp %>% filter(Chrom1 != "Chr14", Chrom2 != "Chr14")

p <- ggplot()+
  geom_rect(data = df.pCp, aes(xmin=ID1-0.5, xmax=ID1+0.5, ymin=ID2-0.5, ymax=ID2+0.5, fill=log2(ReadNum)))+
  geom_rect(data = df.pCp, aes(xmin=ID2-0.5, xmax=ID2+0.5, ymin=ID1-0.5, ymax=ID1+0.5, fill=log2(ReadNum)))+
  scale_fill_gradientn(colours = c("grey80", "#FFA500", "red","darkred"))+
  scale_x_continuous(expand = c(0,0), breaks = df.id$Pos, labels = df.id$Chrom)+
  scale_y_continuous(expand = c(0,0), breaks = df.id$Pos, labels = df.id$Chrom)+
  theme_bw()+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(f_outpdf, p, width = 4, height = 2.5)
write_tsv(df.pCp, f_outtxt)
