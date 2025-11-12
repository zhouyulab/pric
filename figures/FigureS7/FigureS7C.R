rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_EP = "./pRIC/results/TE/EP.pair.TE.txt"
f_outtxt = "./figures/results/FigureS7/FigureS7C/FigureS7C.txt"
f_outpdf = "./figures/results/FigureS7/FigureS7C/FigureS7C.pdf"

df.EP = read_delim(f_EP, "\t")
df.EP.count = df.EP %>% group_by(TEGroup, Group) %>% summarise(num=n()) %>% transform(total = nrow(df.EP)) %>% transform(percent = round(num*100/total, 3))
df.EP.count = df.EP.count %>% transform(label = paste(TEGroup, "(", percent, "%)", sep = ""))
df.EP.count$TEGroup = factor(df.EP.count$TEGroup, levels = c("EP-pair", "enhancer", "promoter", "non-TE"))
p <- ggplot(data = df.EP.count, aes(x=factor(1), y=percent, fill=label)) +
  geom_col(colour = "white") +
  coord_polar(theta = "y", start = 1.65) +
  scale_fill_manual(values = c('#75B7D1FF', '#FF7676FF','#F9D662FF',  '#7CAB7DFF')) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
p

ggsave(f_outpdf, p, width = 4, height = 2.5)
write_tsv(df.EP.count, f_outt)
