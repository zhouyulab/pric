rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_GWAS = "./GWAS/results/GWAS.chimeric.txt"
f_random = "./GWAS/results/GWAS.random.txt"
f_outtxt = "./figures/results/Figure6/Figure6C/Figure6C.txt"
f_outpdf = "./figures/results/Figure6/Figure6C/Figure6C.pdf"

load_df <- function(f_in){
  df = read_delim(f_in, "\t")
  total_reads = df$Num[df$Group=="total"]
  df = df[-1,]
  df = df %>% transform(percent = Num*100/total_reads)
  return(df)
}

df.GWAS = load_df(f_GWAS)
df.random = load_df(f_random) 
df = rbind(df.GWAS, df.random)
df.sum = df %>% group_by(Group) %>% summarise(mean=mean(percent), sd = sd(percent))
df = inner_join(df, df.sum)
df = df %>% transform(zcore = (percent-mean)/sd)

p <- ggplot(data = df, aes(x=Index, y=percent, color=Group))+
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  # geom_line()+
  scale_x_continuous(breaks = c(0, 20, 39), labels = c("-1kb", "junction-site", "1kb"))+
  scale_color_manual(values = c("#EE2C2C", "grey80"))+
  labs(x="")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(colour = "black"))
p
ggsave(f_outpdf, p, width = 4, height = 2.3)
write_tsv(df, f_outtxt)

