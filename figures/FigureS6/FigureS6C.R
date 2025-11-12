rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.chimeric.txt"
f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.EP.chimeric.txt"
f_outtxt = "./figures/results/FigureS6/FigureS6C/FigureS6C.txt"
f_outpdf = "./figures/results/FigureS6/FigureS6C/FigureS6C.pdf"

df.AD = read_delim(f_AD, "\t") %>% filter(Group=="E-P", DonorChrom==AcceptorChrom) %>% transform(Sample="Tetraploid")
df.AA = read_delim(f_AA, "\t") %>% filter(Group=="E-P", DonorChrom==AcceptorChrom) %>% transform(Sample="Diploid")
df.AA = df.AA %>% select(-c("background"))
df = rbind(df.AD, df.AA) %>% transform(len=abs(AcceptorCJS-DonorCJS))
p <- ggplot(data = df, aes(x=log10(len), color=Sample))+
  geom_density(linewidth=0.8)+
  labs(x="distance (log2)") +
  geom_vline(xintercept = c(log10(2800), log10(8000), log10(34000), log10(5000)), color="grey50", linetype="dashed")+
  scale_x_continuous(breaks = c(2, log10(5000), log10(2800), log10(8000), log10(34000), 4, 6, 8), labels = c(2, "5kb","2.8kb", "8kb", "34kb", 4, 6, 8))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(.8, .8),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"))
p              

ggsave(f_outpdf, p, width = 3.5, height = 2.5)
write_tsv(df, f_outtxt)
