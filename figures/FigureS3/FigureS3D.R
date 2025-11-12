rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggpointdensity)

f_rep1 = "./GRO/results/featureCount/PRO_AD_WT_Rep1.FPKM.txt"
f_rep2 = "./GRO/results/featureCount/PRO_AD_WT_Rep2.FPKM.txt"
f_outtxt = "./figures/results/FigureS3/FigureS3D/Figure3D.AD.txt"
f_outpdf = "./figures/results/FigureS3/FigureS3D/Figure3D.AD.pdf"

df.rep1 = read_delim(f_rep1, "\t") %>% select(c("Geneid", "FPKM")) %>% filter(FPKM>0) %>%rename("Rep1"="FPKM") 
df.rep2 = read_delim(f_rep2, "\t") %>% select(c("Geneid", "FPKM")) %>% filter(FPKM>0) %>% rename("Rep2"="FPKM")
df = inner_join(df.rep1, df.rep2)
cor_text = cor(log2(df$Rep1), log2(df$Rep2))
p <- ggplot()+
  geom_pointdensity(data = df, aes(x=log2(Rep1), y=log2(Rep2)), size=0.8)+
  scale_color_gradient(low = "grey80", high = 'red')+
  labs(x="Transcript level in PRO-cap\nRep1 in tetraploid(log2FPKM)", y="Transcript level in PRO-cap\nRep2 in tetraploid(log2FPKM)")+
  annotate("text", x=-2, y=8, label=sprintf("cor = %s", round(cor_text, 2)))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p
ggsave(f_outpdf, p, width = 4, height = 4)
write_tsv(df, f_outtxt)

rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggpointdensity)

f_rep1 = "./GRO/results/featureCount/PRO_AA_WT_Rep1.FPKM.txt"
f_rep2 = "./GRO/results/featureCount/PRO_AA_WT_Rep2.FPKM.txt"
f_outtxt = "./figures/results/FigureS3/FigureS3D/Figure3D.AA.txt"
f_outpdf = "./figures/results/FigureS3/FigureS3D/Figure3D.AA.pdf"

df.rep1 = read_delim(f_rep1, "\t") %>% select(c("Geneid", "FPKM")) %>% filter(FPKM>0) %>%rename("Rep1"="FPKM") 
df.rep2 = read_delim(f_rep2, "\t") %>% select(c("Geneid", "FPKM")) %>% filter(FPKM>0) %>% rename("Rep2"="FPKM")
df = inner_join(df.rep1, df.rep2)
cor_text = cor(log2(df$Rep1), log2(df$Rep2))
p <- ggplot()+
  geom_pointdensity(data = df, aes(x=log2(Rep1), y=log2(Rep2)), size=0.8)+
  scale_color_gradient(low = "grey80", high = 'red')+
  labs(x="Transcript level in PRO-cap\nRep1 in diploid(log2FPKM)", y="Transcript level in PRO-cap\nRep2 in diploid(log2FPKM)")+
  annotate("text", x=-2, y=8, label=sprintf("cor = %s", round(cor_text, 2)))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p
ggsave(f_outpdf, p, width = 4, height = 4)
write_tsv(df, f_outtxt)
