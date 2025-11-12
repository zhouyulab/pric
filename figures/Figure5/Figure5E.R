rm(list = ls())
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(VennDiagram)

f_outtxt = "./figures/results/Figure5/Figure5E/Figure5E.txt"
f_outpdf = "./figures/results/Figure5/Figure5E/Figure5E.pdf"
f_wt = "./pRIC/results/EP/AD_WT_merge/AD_WT_merge.EP.fisher.txt"
f_fl = "./pRIC/results/EP/AD_fl_merge/AD_fl_merge.EP.fisher.txt"

df.wt = read_delim(f_wt, "\t") %>% transform(Sample = "WT") %>% filter(Group=="E-P")
df.fl = read_delim(f_fl, "\t") %>% transform(Sample = "fl")  %>% filter(Group=="E-P")
df.wt = df.wt %>% transform(EP=paste(Enhancer, Promoter, sep = "&"))
df.fl = df.fl %>% transform(EP=paste(Enhancer, Promoter, sep = "&")) 
common_EP = length(unique(c(df.wt$EP, df.fl$EP)))
df.wt = df.wt %>% filter(Reads>3)
df.fl = df.fl %>% filter(Reads>3)

wt_EP = unique(df.wt$EP)
fl_EP = unique(df.fl$EP)
p <-  venn.diagram(
  list(wt_EP, fl_EP),
  category.names = c("WT" , "fl"),
  filename = NULL, cex = 2, cat.cex=1.5, cat.default.pos='text',
  fill = c("#2ec4b6", "#ff9f1c"), cat.pos = c(1,1), col = "transparent", lwd = 2
)

pdf(f_outpdf)
grid.draw(p)
dev.off()


inter = length(intersect(wt_EP, fl_EP))
a = length(wt_EP)
b = length(fl_EP)
phyper(inter-1, b, common_EP-b, a, lower.tail = T)
