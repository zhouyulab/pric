rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

f_EP = "./AA_ChIP/results/conserved/EP.conserved.pair.txt"
f_enhancer = "./AA_ChIP/results/conserved/enhancer.pair.txt"
f_outtxt = "./figures/results/FigureS6/FigureS6F/FigureS6F.txt"
f_outpdf = "./figures/results/FigureS6/FigureS6F/FigureS6F.pdf"

df.EP = read_delim(f_EP, "\t")
conseved_enhancers = unique(df.EP$AD.enhancer[df.EP$Group=="conserved"])

df.enhancer = read_delim(f_enhancer, "\t", col_names = FALSE)
colnames(df.enhancer) <- c("AA.enhancer", "AA.Pos", "AD.enhancer", "AD.Pos", "Length", "identity", "A.ratio", "AD.Ratio")

f_anno = "./AD_ChIP/results/feature/enhancer.group.anno.bed"
df.anno = read_delim(f_anno, "\t", col_names = FALSE) %>% select(c("X4", "X8"))
colnames(df.anno) <- c("AD.enhancer", "Pos")

df.EP = left_join(df.EP, df.anno)
df = df.enhancer %>% filter(AD.enhancer %in% df.EP$AD.enhancer[df.EP$Group=="conserved"]) %>% arrange(identity)
df = left_join(df, df.anno)
df = df[!duplicated(df[, "AD.enhancer"]),]
df$AD.enhancer = factor(df$AD.enhancer, levels = df$AD.enhancer)
p <- ggplot(data = df, aes(x=AD.enhancer, y=identity, size=Length))+
  geom_segment(aes(x = AD.enhancer, xend = AD.enhancer, y = 90, yend = identity),size = .75, show.legend = FALSE, linetype="dashed", color="grey80") +
  geom_point(color="#007BC3")+
  # coord_flip()+
  labs(x="enhancer", y="alignment identity")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black")
  )
p

ggsave(f_outpdf, p, width = 10.5, height = 2)
df = df %>% dplyr::select(c("AA.enhancer", "AD.enhancer", "Length", "identity"))
write_tsv(df, f_outtxt)

