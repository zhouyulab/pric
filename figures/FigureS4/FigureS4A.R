rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

numPts_below_line = function(myVector,slope,x){
  yPt = myVector[x]; b = yPt-(slope*x)
  xPts = 1:length(myVector)
  return(sum(myVector <= (xPts*slope+b)))
}

fopt = function(v) {
  v = sort(v); v[v<0] = 0
  slope = (max(v)-min(v))/length(v)
  x = floor(optimize(numPts_below_line, lower = 1, upper = length(v), myVector=v, slope=slope)$minimum)
  # length(v) - x + 1
  v[x]
}
f_in = "./AA_ChIP/results/fature/enhancer.count.bed"
f_outtxt = "./figures/results/FigureS4/FigureS4A/Figure4A.txt"
f_outpdf = "./figures/results/FigureS4/FigureS4A/Figure4A.pdf"

df <- read_delim(f_in, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
df$H3K27ac = scale(df$X8)
df$DNase_Rep1 = scale(df$X9)
df$DNase_Rep2 = scale(df$X10)

df = df %>% filter(X8>100)
df$MeanCount = rowMeans(df[,c("H3K27ac", "DNase_Rep1", "DNase_Rep2")])
cutoff = fopt(df$MeanCount)

df = df %>% transform(Enhnacer=ifelse(MeanCount>cutoff, "SE", "E"))
df_TE = df %>% filter(Enhnacer == "E") %>% dplyr::arrange(X1,X2)
df_TE$Name = paste("E", 1:nrow(df_TE), sep = "")

df_SE = df %>% filter(Enhnacer == "SE") %>% dplyr::arrange(X1,X2)
df_SE$Name = paste("SE", 1:nrow(df_SE), sep = "")

outdf = rbind(df_TE, df_SE) %>% dplyr::select(c("X1", "X2", "X3", "Name","X4","X6", "X7"))
write_tsv(outdf, f_outtxt, col_names = FALSE)

df.plot = data.frame(values=sort(df$MeanCount), index = 1:length(df$MeanCount))
df.plot = df.plot %>% transform(Group = ifelse(values>cutoff, "SE", "E"))
df.count = df.plot %>% group_by(Group) %>% summarise(num=n()) %>% transform(label=sprintf("%s:%s",Group, num))
df.plot = inner_join(df.plot, df.count) 
p <- ggplot(data = df.plot, aes(x=index, y=values, color=label))+
  geom_hline(yintercept = cutoff, linetype="dashed")+
  geom_point(size=0.8)+
  scale_color_manual(values = c("grey70", "red"))+
  labs(x="", y="mean z-score")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.2, .7),
        axis.text = element_text(colour = "black"))
p
ggsave(f_outpdf,p, width = 2.8, height = 2)
