rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS10/FigureS10C/FigureS10C.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10C/FigureS10C.pdf"

f_anno = "./GWAS/results/GWAS.pos.anno.txt"
df.anno = read_delim(f_anno, "\t")
df.count = df.anno %>% group_by(Pos) %>% summarise(num=n()) %>% transform(Total=nrow(df.anno)) %>%
  transform(Percent = round(num*100/Total, 4))

df.count = df.count %>% transform(label = paste(Pos, "(", Percent, "%)", sep = ""))
p <- ggplot(data = df.count, aes(x=factor(1), y=Percent, fill=label)) +
  geom_col(colour = "white") +
  coord_polar(theta = "y", start = 1.65) +
  scale_fill_manual(values = c('#F9D662FF', '#FF7676FF', '#75B7D1FF', '#7CAB7DFF')) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.title = element_blank()
  )
p

ggsave(f_outpdf, p, width = 4, height = 2.5)
outdf = df.anno %>% dplyr::select(c("Chrom", "Start", "End", "Ref", "Alt", "Reference", "Pos"))
write_tsv(outdf, f_outtxt)
