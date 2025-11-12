rm(list = ls())
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

f_outtxt = "./figures/results/Figure3/Figure3D/Figure3D.txt"
f_outpdf = "./figures/results/Figure3/Figure3D/Figure3D.pdf"

f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.all.fisher.txt"
df.AD = read_delim(f_AD, "\t")
df.AD.count = df.AD %>% group_by(Group) %>% summarise(num=n()) %>% transform(total=nrow(df.AD)) %>% transform(label=paste0(Group, "(", num, ")"))

p_AD <- ggplot(data = df.AD.count, aes(x=factor(1), y=num, fill=label)) +
  geom_col(colour = "white") +
  coord_polar(theta = "y", start = 1.65) +
  scale_fill_brewer(palette = "Accent")+
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
p_AD


f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.all.fisher.txt"
df.AA = read_delim(f_AA, "\t")
df.AA.count = df.AA %>% group_by(Group) %>% summarise(num=n()) %>% transform(total=nrow(df.AA)) %>% transform(label=paste0(Group, "(", num, ")"))

p_AA <- ggplot(data = df.AA.count, aes(x=factor(1), y=num, fill=label)) +
  geom_col(colour = "white") +
  coord_polar(theta = "y", start = 1.65) +
  scale_fill_brewer(palette = "Accent")+
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
p_AA

p <- plot_grid(p_AA, p_AD, labels=c("Diploid", "Tetraploid"))
df = rbind(df.AA.count %>% transform(Group="Diploid"), df.AD.count %>% transform(Group="Tetraploid"))
ggsave(f_outpdf, p, width = 8, height = 2.5)
write_tsv(df, f_outtxt)
