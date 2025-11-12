rm(list = ls())
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(argparse)
library(ggsignif)

science_point <- function(value1, value2, n, label){
  pvalue = ks.test(value1, value2)$p.value
  temp <- as.character(pvalue)
  if("e" %in% base::strsplit(temp, split = "")[[1]]){
    temp <- base::strsplit(temp, split = "e")[[1]]
    num <- round(as.numeric(temp[1]), n)
    index <- as.numeric(temp[2])
    p <- num * 10^index
  }else{
    p <- round(pvalue, n)
  }
  Plable = sprintf("P(%s) = %s", label, p)
  return(Plable)
}

load_EP_gene <- function(f_in, cutoff=0){
  EP_gene = c()
  df = read_delim(f_in, "\t") %>% filter(Group=="E-P") %>% filter(Reads>cutoff)
  for(i in 1:nrow(df)){
    gene = gsub("pro", "", df$Promoter[i])
    EP_gene = c(EP_gene, gene)
  }
  return(EP_gene)
}

plot_EP <- function(df_RNA, gene_li, samp){
  df_RNA = df_RNA %>% transform(Type=ifelse(GeneName %in% gene_li, "T", "F")) %>% filter(FPKM>1)
  df_RMA_count = df_RNA %>% group_by(Type) %>% summarise(n=n())
  df_RNA = left_join(df_RNA, df_RMA_count) %>% transform(label = paste(Type, "(", n, ")", sep = ""))
  
  
  T_fpkm = df_RNA$FPKM[df_RNA$Type=="T"]
  F_fpkm = df_RNA$FPKM[df_RNA$Type=="F"]
  pvalue = science_point(log2(T_fpkm), log2(F_fpkm), 3, "T vs F")
  p <- ggplot(data = df_RNA, aes(x=log2(FPKM), color=label))+
    stat_ecdf(linewidth=0.8)+
    scale_color_brewer(palette = "Set1")+
    annotate("text", x=5, y=0.5, label=pvalue)+
    scale_x_continuous(limits = c(0, 7))+
    labs(y="Cumulative fraction", x=sprintf("%s gene expression (log2)", samp))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "right",
          axis.text = element_text(colour = "black"))
  
  return(p)
}

f_fl = "./pRIC/results/EP/AD_fl_merge/AD_fl_merge.EP.fisher.txt"
f_RNA = "./AD_RNA/results/featureCount/RNA.FPKM.txt"
f_outtxt = "./figures/results/FigureS8/FigureS8A/FigureS8A.txt"
f_outpdf = "./figures/results/FigureS8/FigureS8A/FigureS8A.pdf"

df_RNA = read_delim(f_RNA, "\t")
df_RNA = df_RNA %>% transform(WT = (WT_Rep1+WT_Rep2)/2, fl = (fl_Rep1+fl_Rep2)/2)
fl_gene = load_EP_gene(f_fl)

df.RNA.fl = df_RNA %>% dplyr::select(c("GeneName", "fl")) %>% dplyr::rename("FPKM"="fl")
p_fl = plot_EP(df.RNA.fl, fl_gene, "fl")

ggsave(f_outpdf, p_fl, width = 3.5, height = 2)
df.RNA.fl = df.RNA.fl %>% transform(Type=ifelse(GeneName %in% fl_gene, "T", "F")) %>% filter(FPKM>1)
write_tsv(df.RNA.fl,f_outtxt)
