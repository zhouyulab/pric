rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(argparse)

f_outtxt = "./figures/results/Figure3/Figure3E/Figure3E.txt"
f_outpdf = "./figures/results/Figure3/Figure3E/Figure3E.pdf"

f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.nearest.txt"
df.AD = read_delim(f_AD, "\t") %>% transform(Sample="Tetraploid")
df.AD.count = df.AD %>% group_by(Position) %>% summarise(n=n()) %>% transform(percent=round(n/nrow(df.AD), 2))
p_AD <- ggplot(data = df.AD.count, aes(x=Position, y=percent, fill=Position))+
  geom_bar(stat = "identity", width=0.6)+
  geom_text(data = df.AD.count, aes(x=Position, y=percent+0.1, label=paste(percent*100, "%", sep = "")))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y="percent")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"))
p_AD


f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.EP.nearest.txt"
df.AA = read_delim(f_AA, "\t") %>% transform(Sample="Diploid")
df.AA.count = df.AA %>% group_by(Position) %>% summarise(n=n()) %>% transform(percent=round(n/nrow(df.AA), 2))
p_AA <- ggplot(data = df.AA.count, aes(x=Position, y=percent, fill=Position))+
  geom_bar(stat = "identity", width=0.6)+
  geom_text(data = df.AA.count, aes(x=Position, y=percent+0.1, label=paste(percent*100, "%", sep = "")))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y="percent")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"))
p_AA

p <- plot_grid(p_AA, p_AD, labels=c("Diploid", "Tetraploid"))
df = rbind(df.AA.count %>% transform(Group="Diploid"), df.AD.count %>% transform(Group="Tetraploid"))
ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
