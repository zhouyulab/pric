rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/FigureS6/FigureS6A/FigureS6A.txt"
f_outpdf = "./figures/results/FigureS6/FigureS6A/FigureS6A.pdf"

f_AD_promoter = "./AD_ChIP/results/feature/promoter.bed"
f_AD_enhancer = "./AD_ChIP/results/feature/enhancer.group.bed"
df.AD.promoter = read_delim(f_AD_promoter, "\t", col_names = FALSE) %>% transform(length = X3-X2, feature="Promoter", Speceis="Tetraploid")
df.AD.promoter = df.AD.promoter %>% select(c("length", "feature", "Speceis"))
df.AD.enhancer = read_delim(f_AD_enhancer, "\t", col_names = FALSE) %>% transform(length = X3-X2, feature="Enhancer", Speceis="Tetraploid")
df.AD.enhancer = df.AD.enhancer %>% select(c("length", "feature", "Speceis"))

f_AA_promoter = "./AA_ChIP/results/feature/promoter.bed"
f_AA_enhancer = "./AA_ChIP/results/feature/enhancer.group.bed"
df.AA.promoter = read_delim(f_AA_promoter, "\t", col_names = FALSE) %>% transform(length = X3-X2, feature="Promoter", Speceis="Diploid")
df.AA.promoter = df.AA.promoter %>% select(c("length", "feature", "Speceis"))
df.AA.enhancer = read_delim(f_AA_enhancer, "\t", col_names = FALSE) %>% transform(length = X3-X2, feature="Enhancer", Speceis="Diploid")
df.AA.enhancer = df.AA.enhancer %>% select(c("length", "feature", "Speceis"))

df = Reduce(rbind, list(df.AD.promoter, df.AD.enhancer, df.AA.promoter, df.AA.enhancer))
df.count = df %>% group_by(feature, Speceis) %>% summarise(num=n())
df.info = df %>% group_by(feature, Speceis) %>% summarise(meanLen=mean(length), medianLen = median(length),
                                                          minLen=min(length), maxLen = max(length))
p <- ggplot(data = df, aes(x=Speceis, y=log2(length), fill=Speceis))+
  facet_wrap(vars(feature), scales = "free")+
  geom_violin()+
  scale_fill_manual(values = c("#CB87B4", "#088BBE"))+
  # geom_text(data = df.count, aes(x=Speceis, y=c(8.5, 7.5), label=num))+
  geom_boxplot(width=0.3, outlier.shape = NA)+
  labs(x="", y="length (log2)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"))
p

ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
