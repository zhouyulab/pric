rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS10/FigureS10A/FigureS10A.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10A/FigureS10A.pdf"
f_indir= "./GWAS/data"

f_2023_NG = file.path(f_indir, "2023_NG_eQTL/results/2023_NG_GWAS_info.txt")
df.2023.NG = read_delim(f_2023_NG, "\t")
df.2023.NG = df.2023.NG[!duplicated(df.2023.NG[, c("Chrom", "Start", "End", "Ref", "Alt")]),]
df.eQTL = df.2023.NG %>% group_by(Type) %>% summarise(num=n()) %>% transform(total=nrow(df.2023.NG))
df.eQTL$Type = factor(df.eQTL$Type, labels = c("cis", "trans"))
df.eQTL = df.eQTL %>% transform(percent = round(num*100/total, 4)) %>% transform(label=paste(percent, "%", sep = "")) %>% arrange(percent)
df.eQTL = df.eQTL %>% transform(pos=cumsum(rev(percent)))
p <- ggplot(data = df.eQTL, aes(x=factor(1), y=percent, fill=Type))+
  geom_bar(stat = "identity", color="white")+
  geom_text(aes(x=factor(1), y=rev(pos), label=label))+
  scale_fill_manual(values = c("#008856", "#9999FF"))+
  coord_polar("y", start = 0)+
  theme_void()
p

ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df.eQTL[, c(1:4)], f_outtxt)
