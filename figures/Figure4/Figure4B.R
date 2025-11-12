rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_fl = "./pRIC/results/EP/AD_fl_merge/EP_flank_MFE.txt"
f_WT = "./pRIC/results/EP/AD_WT_merge/EP_flank_MFE.txt"
f_AA = "./pRIC/results/EP/A_WT_merge/EP_flank_MFE.txt"
f_outtxt = "./figures/results/Figure4/Figure4B/Figure4B.txt"
f_outpdf = "./figures/results/Figure4/Figure4B/Figure4B.pdf"

df.fl = read_delim(f_fl, "\t") %>% transform(Sample="fl")
df.WT = read_delim(f_WT, "\t") %>% transform(Sample="WT")
df.AA = read_delim(f_AA, "\t") %>% transform(Sample="AA")

df = Reduce(rbind, list(df.fl, df.WT, df.AA))
df$Sample = factor(df$Sample, levels = c("AA", "WT", "fl"))
p <- ggplot(data = df, aes(x=Sample, y=MFE, fill=Group))+
  stat_boxplot(geom = "errorbar", width=0.3, position=position_dodge(0.7)) +
  geom_boxplot(outlier.shape = NA, width=0.6, position=position_dodge(0.7))+
  scale_y_continuous(limits = c(-50, 0))+
  scale_fill_brewer(palette = "Pastel2")+
  labs(x="", y="MFE")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "right")
p
ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
