rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS10/FigureS10D/FigureS10D.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10D/FigureS10D.pdf"

flanks = c(200, 500, 1000, 2000, 5000, 10000)
flank_li = list()
for(flank in flanks){
  f_in = sprintf("./GWAS/results/eQTL/anno_merge_eQTL_chimeric_%s.txt", flank)
  df = read_delim(f_in, "\t", col_names = FALSE)
  colnames(df) =  c("Chrom", "Start", "End", "Type", "GeneName", "Source", "Table", "RefID","ReadNum", "ingene", "pos")
  df = df %>% filter(ReadNum>0) %>% group_by(Type) %>% summarise(n=n())
  df = df %>% transform(Flank=flank)
  flank_li[[flank]] = df
}
dat = Reduce(rbind, flank_li)
dat$Flank = factor(dat$Flank, levels = flanks, labels = c("0.2kb", "0.5kb", "1kb", "2kb", "5kb", "10kb"))
p <- ggplot(data = dat, aes(x=Flank, y=n, fill=Type))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(Type), scales = "free")+
  scale_fill_manual(values = c("#7CCA89FF", "#508432FF"))+
  labs(x="", y="epistasis count") +
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black"))
p
ggsave(f_outpdf, p, width = 5, height = 2.6)

write_tsv(dat, f_outtxt)
