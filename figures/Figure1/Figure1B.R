rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_wt = "./pRIC/results/QC/pCp/AD_WT_merge.plot.txt"
f_fl = "./pRIC/results/QC/pCp/AD_fl_Rep1.plot.txt"
f_outtxt = "./figures/results/Figure1/Figure1B/Figure1B.txt"
f_outpdf = "./figures/results/Figure1/Figure1B/Figure1B.pdf"

df_wt <- read_delim(f_wt, delim = "\t", escape_double = FALSE,  trim_ws = TRUE) %>% transform(Sample="WT")
df_fl <- read_delim(f_fl, delim = "\t", escape_double = FALSE,  trim_ws = TRUE) %>% transform(Sample="fl")
df = Reduce(rbind, list(df_wt, df_fl))
df$Sample = factor(df$Sample, levels=c("WT", "fl"))
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
