rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggsignif)

f_TF = "./plantTFDB/results/GWAS/TF_contain_GWAS.txt"
f_random = "./plantTFDB/results/GWAS/TF_contain_GWAS.random.txt"
f_outtxt = "./figures/results/FigureS10/FigureS10H/FigureS10H.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10H/FigureS10H.pdf"


df.TF = read_delim(f_TF, "\t")
df.random = read_delim(f_random, "\t")
df = rbind(df.TF, df.random)
df = df %>% filter(GWASnum>0)
df.count = df %>% group_by(Group) %>% summarise(num=n())

p <- ggplot(data = df, aes(x=Group, y=ratio))+
  geom_boxplot(aes(fill=Group), width=0.6, outlier.shape = NA, notch=TRUE, linetype="dashed", position=position_dodge(0.7))+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill=Group), width=0.6, outlier.shape = NA, notch=TRUE, position=position_dodge(0.7))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,color="black")+
  geom_text(data = df.count, aes(x=Group, y=-4, label=num))+
  scale_y_continuous(limits = c(-4, 60))+
  geom_signif(comparisons = list(c("random", "TF")), 
              test = t.test, y_position = c(50),map_signif_level = T, 
              color = "black",textsize = 4,tip_length=0.005) + 
  scale_fill_manual(values = c("#5773CC", "#FFB900"))+
  labs(x="")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"))

p
ggsave(f_outpdf, p, width = 3, height = 4)
write_tsv(df, f_outtxt)
