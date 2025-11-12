rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS10/FigureS10B/FigureS10B.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10B/FigureS10B.pdf"
f_indir= "./GWAS/data"

f_2017_NG = file.path(f_indir, "2017_NG/results/2017_NG_GWAS_info.txt")
df.2017.NG = read_delim(f_2017_NG, "\t")
df.2017.NG = df.2017.NG[!duplicated(df.2017.NG[, c("Chrom", "Start", "End", "Ref", "Alt")]),]

f_2018_NG = file.path(f_indir, "2018_NG/results/2018_NG_GWAS_info.txt")
df.2018.NG = read_delim(f_2018_NG, "\t")
df.2018.NG = df.2018.NG[!duplicated(df.2018.NG[, c("Chrom", "Start", "End", "Ref", "Alt")]),]

f_2021_PBJ = file.path(f_indir, "2021_PBJ/results/2021_PBJ_GWAS_info.txt")
df.2021.PBJ = read_delim(f_2021_PBJ, "\t")
df.2021.PBJ = df.2021.PBJ[!duplicated(df.2021.PBJ[, c("Chrom", "Start", "End", "Ref", "Alt")]),]

f_2023_JGG = file.path(f_indir, "2023_JGG_breed/results/2023_JGG_GWAS_info.txt")
df.2023.JGG = read_delim(f_2023_JGG, "\t")
df.2023.JGG = df.2023.JGG[!duplicated(df.2023.JGG[, c("Chrom", "Start", "End", "Ref", "Alt")]),]

f_2023_NAR = file.path(f_indir, "2023_NAR_cottonMD/results/2023_NAR_GWAS_info.txt")
df.2023.NAR = read_delim(f_2023_NAR, "\t")
df.2023.NAR = df.2023.NAR[!duplicated(df.2023.NAR[, c("Chrom", "Start", "End", "Ref", "Alt")]),]

df = Reduce(rbind, list(df.2017.NG, df.2018.NG, df.2021.PBJ, df.2023.JGG, df.2023.NAR))
df.GWAS = df %>% group_by(Reference) %>% summarise(num=n()) %>% transform(total=nrow(df))
df.GWAS = df.GWAS %>% transform(percent = round(num*100/total, 4)) %>% transform(label=paste(percent, "%", sep = "")) %>% arrange(percent)
df.GWAS$Reference = factor(df.GWAS$Reference, levels = c("2017_NG", "2018_NG", "2021_PBJ", "2023_NAR", "breed"),
                           labels = c("2017_NG", "2018_NG", "2021_PBJ", "2023_NAR", "2023_JGG"))
df.GWAS = df.GWAS %>% arrange(Reference)
p <- ggplot(data = df.GWAS, aes(x=factor(1), y=percent, fill=Reference))+
  geom_bar(stat = "identity", color="white")+
  # geom_text(aes(x=factor(1), y=rev(pos), label=label))+
  scale_fill_manual(values = c("#E62932", "#0099FF", "#008856", "#9999FF", "#FF794B"))+
  coord_polar("y", start = 0)+
  theme_void()
p

ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df.GWAS, f_outtxt)
