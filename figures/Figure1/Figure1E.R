rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_in = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge.ChimericStatis.txt"
f_outtxt = "./figures/results/Figure1/Figure1E/Figure1E.txt"
f_outpdf = "./figures/results/Figure1/Figure1E/Figure1E.pdf"

df = read_delim(f_in, "\t", col_names = FALSE)
colnames(df) <- c("Group", "ReadNum", "Sample")
df$Total = sum(df$ReadNum)
df = df %>% transform(Percent=round(ReadNum*100/Total, 3)) %>% transform(label=paste(Percent, "%", sep = ""))
df = df %>% arrange(Percent) %>% transform(Pos=cumsum(Percent)-0.5*Percent)
df$Group = factor(df$Group, levels = c("None-None", "intra-gene", "gene-None","inter-gene"),
                  labels = c("IGR-IGR", "gene-intra","gene-intIGR", "gene-inter"))
df = df %>% transform(label2 = paste(Group, "(", label, ")", sep = "")) %>% arrange(Group)
df$label2 = factor(df$label2)
p <- ggplot(data = df, aes(x=factor(1), y=Percent, color=label2, fill=label2)) +
  geom_col(colour = "white") +
  # geom_text(aes(x=factor(1), y=Pos, label=label), color="black")+
  coord_polar(theta = "y", start = 1.65) +
  scale_fill_brewer(palette = "Set2")+
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
p
ggsave(f_outpdf, p, width = 4, height = 2)
df = df[,c(1:5)]
write_tsv(df, f_outtxt)
