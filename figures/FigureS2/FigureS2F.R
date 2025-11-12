rm(list=ls())
library(readr)
library(dplyr)
library(ggplot2)
annotate_pvalue <- function(pvalue){
  if(pvalue < 0.001){
    annotation <- "***"
  }else{
    if(pvalue < 0.01){
      annotation <- "**"
    }else{
      if(pvalue < 0.05){
        annotation <- "*"
      }else{
        annotation <- ""
      }
    }
  }
  return(annotation)
}


f_DEseq = "./AD_RNA/results/featureCount/RNA.DEseq.res"
f_wt_chimeric = "./pRIC/results/analysis/AD_WT_merge/AD_WT_merge.intra.Chimeric.txt"
f_fl_chimeric = "./pRIC/results/analysis/AD_fl_merge/AD_fl_merge.intra.Chimeric.txt"
f_outtxt = "./figures/results/FigureS2/FigureS2F/FigureS2F.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2F/FigureS2F.pdf"

df.DEseq = read_delim(f_DEseq, "\t")
df.wt.chimeric = read_delim(f_wt_chimeric, "\t") %>% transform(Sample="WT")
df.wt.chimeric = df.wt.chimeric %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% filter(GeneType=="protein") 
df.fl.chimeric = read_delim(f_fl_chimeric, "\t") %>% transform(Sample="fl")
df.fl.chimeric = df.fl.chimeric %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = "[|]") %>% filter(GeneType=="protein")

df1 = df.DEseq %>% inner_join(df.wt.chimeric)
df2 = df.DEseq %>% inner_join(df.fl.chimeric)
df = rbind(df1, df2)

df.pvalue = data.frame()
for(change in unique(df$change)){
  wt_value = df$ReadNum[df$change==change & df$Sample=="WT"]
  fl_value = df$ReadNum[df$change==change & df$Sample=="fl"]
  pvalue = t.test(log2(wt_value), log2(fl_value))$p.value
  df.pvalue = rbind(df.pvalue, data.frame(chang=change, pvalue=annotate_pvalue(pvalue)))
}
df.pvalue$Sample="WT"
df$Sample = factor(df$Sample, levels = c("WT", "fl"))
p <- ggplot(data = df, aes(x=change, y=log2(ReadNum), fill=Sample))+
  geom_boxplot(aes(fill=Sample), width=0.6, outlier.shape = NA, notch=TRUE, linetype="dashed", position=position_dodge(0.7))+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill=Sample), width=0.6, outlier.shape = NA, notch=TRUE, position=position_dodge(0.7))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,color="black", position=position_dodge(0.7))+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,color="black", position=position_dodge(0.7))+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(limits = c(0, 12))+
  labs(x="", y="chimeric reads (log2)")+
  geom_text(data = df.pvalue, aes(x=chang, y=11, label=pvalue))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"))
p
ggsave(f_outpdf, p, width = 4, height = 3)
write_tsv(df, f_outtxt)
