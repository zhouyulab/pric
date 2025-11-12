rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
f_in = "./AD_RNA/results/featureCount/RNA.DEseq.res"
f_outtxt = "./figures/results/FigureS2/FigureS2E/FigureS2E.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2E/FigureS2E.pdf"

df = read_delim(f_in, "\t")
df_count = df %>% group_by(change) %>% dplyr::summarise(n=n()) %>% transform(label = sprintf("%s (%s)", change, n))
p <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(pvalue), color=change))+
  geom_point(size=0.8)+
  scale_color_manual(values = c("#3952a3", "grey80", "#ef2224"), label = df_count$label)+
  labs(x= "expression foldchange (log2)", y="log10(p-value)")+
  scale_x_continuous(limits = c(-15, 15))+
  scale_x_continuous(limits = c(-15, 15))+
  geom_vline(xintercept = c(-1, 1), linetype="dashed",color="grey50")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed",color="grey50")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = c(.2, .8))
p
ggsave(f_outpdf, p, width=4, height =3)
write_tsv(df, f_outtxt)
