rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS8/FigureS8C/FigureS8C.VIGS1.txt"
f_outpdf = "./figures/results/FigureS8/FigureS8C/FigureS8C.VIGS1.pdf"
f_protein = "./VIGS_RNA/results/featureCount/VIGS1.DEseq.txt"
f_lncRNA = "./VIGS_RNA/results/lncRNA/VIGS1.DEseq.txt"

df.protein = read_delim(f_protein, "\t")
df.lncRNA = read_delim(f_lncRNA, "\t")
df = rbind(df.protein, df.lncRNA)

df$change = "nochange"
df$change[df$padj<0.05 & df$log2FoldChange> log2(2)] <- "up"
df$change[df$padj<0.05 & df$log2FoldChange < log2(1/2)] <- "down"

df_count = df %>% group_by(change) %>% dplyr::summarise(n=n()) %>% transform(label = sprintf("%s (%s)", change, n))
p <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), color=change))+
  geom_point(size=0.8)+
  scale_color_manual(values = c("#3952a3", "grey80", "#ef2224"), label = df_count$label)+
  labs(x= "expression foldchange (log2)", y="log10(p-value)")+
  scale_x_continuous(limits = c(-7, 7))+
  # scale_y_continuous(limits = c(0, 200))+
  scale_size_manual(values = c(0.8, 3.5))+
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


rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS8/FigureS8C/FigureS8C.VIGS2.txt"
f_outpdf = "./figures/results/FigureS8/FigureS8C/FigureS8C.VIGS2.pdf"
f_protein = "./VIGS_RNA/results/featureCount/VIGS2.DEseq.txt"
f_lncRNA = "./VIGS_RNA/results/lncRNA/VIGS2.DEseq.txt"

df.protein = read_delim(f_protein, "\t")
df.lncRNA = read_delim(f_lncRNA, "\t")
df = rbind(df.protein, df.lncRNA)

df$change = "nochange"
df$change[df$padj<0.05 & df$log2FoldChange> log2(2)] <- "up"
df$change[df$padj<0.05 & df$log2FoldChange < log2(1/2)] <- "down"

df_count = df %>% group_by(change) %>% dplyr::summarise(n=n()) %>% transform(label = sprintf("%s (%s)", change, n))
p <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), color=change))+
  geom_point(size=0.8)+
  scale_color_manual(values = c("#3952a3", "grey80", "#ef2224"), label = df_count$label)+
  labs(x= "expression foldchange (log2)", y="log10(p-value)")+
  scale_x_continuous(limits = c(-7, 7))+
  # scale_y_continuous(limits = c(0, 200))+
  scale_size_manual(values = c(0.8, 3.5))+
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
