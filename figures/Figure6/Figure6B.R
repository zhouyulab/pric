rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_outtxt = "./figures/results/Figure6/Figure6B/Figure6B.txt"
f_outpdf = "./figures/results/Figure6/Figure6B/Figure6B.pdf"
f_in = "./GWAS/results/eQTL.plot.txt"

df = read_delim(f_in, "\t")
df.sum = df %>% group_by(Group) %>% summarise(mean=mean(Percent), sd = sd(Percent))
df = inner_join(df, df.sum)
df = df %>% transform(zcore = (Percent-mean)/sd)
p <- ggplot(data = df, aes(x=Index, y=Percent, color=Group))+
  geom_smooth(method = "loess", span = 0.4, se = TRUE) +
  scale_x_continuous(breaks = c(0, 20, 39), labels = c("eQTL-site", "1kb", "2kb"))+
  scale_color_manual(values = c("#EE2C2C", "grey70"))+
  labs(x="")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(colour = "black"))
p

ggsave(f_outpdf, p, width = 4, height = 2.3)
write_tsv(df, f_outtxt)
