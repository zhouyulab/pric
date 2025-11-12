rm(list=ls())
library(readr)
library(dplyr)
library(ggplot2)

f_in = "./pRIC/results/QC/pCp/A_WT_merge.plot.txt"
f_outtxt = "./figures/results/FigureS2/FigureS2A/FigureS2A.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2A/FigureS2A.pdf"


df <- read_delim(f_in, delim = "\t", escape_double = FALSE,  trim_ws = TRUE) %>% transform(Sample="WT")
p <- ggplot(data = df, aes(x=indx, y=C, color=Sample))+
  geom_line()+
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20), 
                     labels = c("-20 bp", "-10 bp", "Junction", "10bp", "20 bp"))+
  scale_y_continuous(breaks = c(20, 40, 60, 80, 100))+
  scale_color_brewer(palette = "Set1")+
  labs(x="", y="Cytosine (%)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, color = "black", hjust = 1),
        axis.text.y = element_text(color = "black", hjust = 1),
        legend.key.size = unit(8, "pt"),
        legend.title = element_blank(),
        legend.position = c(.8,.7))
p
ggsave(f_outpdf, p, width = 3,height = 2.5)
write_tsv(df, f_outtxt)
