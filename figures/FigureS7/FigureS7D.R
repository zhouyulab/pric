rm(list = ls())
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
        annotation <- "NS"
      }
    }
  }
  return(annotation)
}
load_random_pair <- function(){
  sample_li = c("WT_Rep1", "WT_Rep2", "fl_Rep1", "fl_Rep2")
  df.random = data.frame()
  for(samp in sample_li){
    f_random_MFE = sprintf("./pRIC/results/TE/chimeric/%s_250_random_EP.MFE.txt", samp)
    df.random.MFE = read_delim(f_random_MFE, "\t") %>% filter(Group=="EP")
    colnames(df.random.MFE) <- c("Promoter","PromoterChrom","PromoterSite","Enhancer","EnhancerChrom", "EnhancerSite", "Group", "MFE")
    f_random_TE = sprintf("./pRIC/results/TE/chimeric/%s_random_EP.chimeric.250.TE.txt", samp)
    df.random.TE = read_delim(f_random_TE, "\t")
    df.random.temp = inner_join(df.random.MFE, df.random.TE) %>% transform(pos="random")
    df.random = rbind(df.random, df.random.temp)
  }
  return(df.random)
}
load_EP_pair <- function(){
  f_MFE = "./pRIC/results/EP/AD_pCp_merge/EP_flank_MFE.txt"
  df.MFE = read_delim(f_MFE, "\t") %>% filter(Group=="EP")
  f_TE = "./pRIC/results/TE/chimeric/EP.chimeric.250.TE.txt"
  df.TE = read_delim(f_TE, "\t")
  df = inner_join(df.MFE, df.TE)
  df = df %>% transform(pos="EP")
  return(df)
}
group_TE = function(x, RNA_TE, DNA_TE, target_TE){
  enhancer_te = x[11]
  te_li = strsplit(enhancer_te, split = "[|]") %>% unlist() %>% unique()
  if(length(te_li) == 1){
    if(te_li %in% target_TE){return(te_li)}
    if(te_li %in% RNA_TE & !te_li %in% target_TE){return("RNA-other")}
    if(te_li %in% DNA_TE & !te_li %in% target_TE){return("DNA-other")}
  }
  if("LTR/Gypsy" %in% te_li){return("LTR/Gypsy")}
  if("LINE/L1" %in% te_li){return("LINE/L1")}
  if("Simple_repeat" %in% te_li){return("Simple_repeat")}
  if("Low_complexity" %in% te_li){return("Low_complexity")}
  if("LTR/Copia" %in% te_li){return("LTR/Copia")}
  if("Unknown" %in% te_li){return("Unknown")}
  if("DNA/MULE-MuDR" %in% te_li){return("DNA/MULE-MuDR")}
  if("DNA/hAT-Ac" %in% te_li){return("DNA/hAT-Ac")}
  if(length(intersect(RNA_TE, te_li))>0){return("RNA-other")}
  if(length(intersect(DNA_TE, te_li))>0){return("DNA-other")}
  return("non-TE")
}
calculate_pvalue <- function(df){
  TE_group_li = unique(df$Group2)
  df.pvalue = data.frame()
  for(te in TE_group_li){
    EP_MFE = df$MFE[df$Group2==te & df$pos=="EP"]
    random_MFE = df$MFE[df$Group2==te & df$pos=="random"]
    pvalue = t.test(EP_MFE, random_MFE)$p.value
    panno = annotate_pvalue(pvalue)
    df.temp = data.frame(Group2=te, pvalue = pvalue, panno=panno)
    df.pvalue = rbind(df.pvalue, df.temp)
  }
  return(df.pvalue)
}

f_outtxt = "./figures/results/FigureS7/FigureS7D/FigureS7D.txt"
f_outpdf = "./figures/results/FigureS7/FigureS7D/FigureS7D.pdf"

df = load_EP_pair()
df.random = load_random_pair()
df = rbind(df, df.random)
RNA_TE = c("Unknown", "LTR/Gypsy","Simple_repeat", "Low_complexity", "LINE/L1", "LTR/Copia", "SINE/tRNA", "tRNA", "rRNA", "LTR/Caulimovirus", "LTR/ERV4", "Satellite")
DNA_TE = c("DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/hAT", "DNA/hAT-Tag1", "DNA/hAT-Ac", "DNA", "DNA/hAT-Tip100", "RC/Helitron", "DNA/CMC-EnSpm", "DNA/TcMar", "DNA/Maverick")
target_TE = c("Unknown", "LTR/Gypsy", "LTR/Copia", "DNA/MULE-MuDR", "DNA/hAT-Ac", "LINE/L1", "Simple_repeat", "Low_complexity")

df$Group2 = apply(df,1,group_TE, RNA_TE, DNA_TE,target_TE)
df = df %>% filter(Group2 != "non-TE")
df$Group2 = factor(df$Group2, levels = c("Simple_repeat", "LTR/Gypsy", "LINE/L1", "LTR/Copia","RNA-other",
                                         "DNA/MULE-MuDR", "DNA/hAT-Ac", "DNA-other", "Low_complexity", "Unknown"))
df.pvalue = calculate_pvalue(df) %>% transform(pos="EP")
p <- ggplot(data = df, aes(x=Group2, y=MFE, fill=pos))+
  stat_boxplot(geom = "errorbar", width=0.2, position=position_dodge(0.7)) +
  geom_boxplot(outlier.shape = NA, width=0.6, position=position_dodge(0.7))+
  scale_y_continuous(limits = c(-50,-0))+
  geom_text(data = df.pvalue, aes(x=Group2, y=-5, label=panno))+
  scale_fill_manual(values = c("#990099","#009900"))+
  labs(x="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top")
p

ggsave(f_outpdf, p, width =5, height = 3.5)
write_tsv(df, f_outtxt)
