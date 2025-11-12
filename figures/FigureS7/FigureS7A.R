rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)

load_streme <- function(f_sequence, f_bed, f_streme){
  df.sequence = read_delim(f_sequence, "\t") %>% na.omit()
  df.bed = read_delim(f_bed, "\t", col_names = FALSE)
  colnames(df.bed) <- c("Chrom", "Start", "End", "Name", "Pos", "Strand")
  df.bed = df.bed %>% transform(ID=paste(Name, "::", Chrom, ":", Start, "-", End, "(", Strand, ")", sep = ""))
  df.sequence.count = df.sequence %>% group_by(motif_ID, `motif_P-value`) %>% summarise(num=n()) %>% transform(total = nrow(df.bed))
  
  streme=file(f_streme, open="r")
  n=1
  motif.list  = list()
  while ( TRUE ) {
    line = readLines(streme, n = 1)
    if(length(line) == 0){
      break
    }
    if(startsWith(line,  "ALPHABET")){
      seq = strsplit(line, split = "= ")[[1]][2]
      seq = strsplit(seq, split = "")[[1]]
    }
    if(startsWith(line, "MOTIF")){
      # motif = strsplit(strsplit(line, split = " ")[[1]][2], " ")[[1]][1]
      motif = strsplit(line, split = " ")[[1]][2]
      info = readLines(streme, n = 1)
      alength = strsplit(str_match_all(info,pattern = "alength=.*w")[[1]], " ")[[1]][2]
      w = strsplit(str_match_all(info,pattern = "w=.*nsites")[[1]], " ")[[1]][2]
      nsites = strsplit(str_match_all(info,pattern = "nsites=.*E")[[1]], " ")[[1]][2]
      E = strsplit(str_match_all(info,pattern = "E=.*")[[1]], "= ")[[1]][2]
      nums = readLines(streme, n = as.numeric(w))
      num = c()
      for(x in nums){
        a = trimws(x, which = "both")
        a = strsplit(a, " ")[[1]]
        num = c(num, a)
      }
      seqmatrix <- t(matrix(as.numeric(num), ncol=4,  byrow = T))
      rownames(seqmatrix) <- seq
      motif.list[[motif]] = list(alength=alength, w=w, nsites=nsites, E=E,  matrix=seqmatrix)
      n = n+1+as.numeric(w)
    }else{
      n = n+1
    }
  }
  close(streme)
  
  df.sequence.count = df.sequence.count %>% transform(Evalue = 0, nsites =0)
  for(name in names(motif.list)){
    nsite = motif.list[[name]]$nsites
    E_value = motif.list[[name]]$E
    df.sequence.count$Evalue[df.sequence.count$motif_ID==name] = E_value
    df.sequence.count$nsites[df.sequence.count$motif_ID==name] = as.numeric(nsite)
  }
  df.sequence.count = df.sequence.count %>% transform(percent = round(nsites*100/total, 2)) %>% arrange(desc(percent)) 
  df.sequence.count$motif_ID = factor(df.sequence.count$motif_ID)
  return(res = list(df=df.sequence.count, motif=motif.list))
}

f_outtxt = "./figures/results/FigureS7/FigureS7A/FigureS7A.txt"
f_outpdf = "./figures/results/FigureS7/FigureS7A/FigureS7A.pdf"

f_enhancer_sequence = "./pRIC/results/EP/AD_pCp_merge/streme_enhancer/sequences.tsv"
f_enhancer_bed = "./pRIC/results/EP/AD_pCp_merge/enhancer.bed"
f_enhancer_streme = "./pRIC/results/EP/AD_pCp_merge/streme_enhancer/streme.txt"
enhancer_res = load_streme(f_enhancer_sequence, f_enhancer_bed, f_enhancer_streme)
df.enhancer = enhancer_res$df[1:20,]
df.enhancer$motif_ID = factor(df.enhancer$motif_ID, levels = df.enhancer$motif_ID)
p_enhancer <- ggplot(data = df.enhancer, aes(x=motif_ID, y=percent))+
  geom_bar(stat = "identity", fill="orange")+
  labs(x="", y="Percent")+
  geom_text(data = df.enhancer, aes(x=motif_ID, y=percent, label=motif_ID))+
  scale_y_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0, 10, 20, 30, 40), labels = c(0, 10, 20, 30, 40))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_blank()
  )
p_enhancer

f_promoter_sequence = "./pRIC/results/EP/AD_pCp_merge/streme_promoter/sequences.tsv"
f_promoter_bed = "./pRIC/results/EP/AD_pCp_merge/promoter.bed"
f_promoter_streme = "./pRIC/results/EP/AD_pCp_merge/streme_promoter/streme.txt"
promoter_res = load_streme(f_promoter_sequence, f_promoter_bed, f_promoter_streme)
df.promoter = promoter_res$df[1:20,]
df.promoter$motif_ID = factor(df.promoter$motif_ID, levels = df.promoter$motif_ID)
p_promoter <- ggplot(data = df.promoter, aes(x=motif_ID, y=-percent))+
  geom_bar(stat = "identity", fill="orange")+
  labs(x="", y="Percent")+
  geom_text(data = df.promoter, aes(x=motif_ID, y=-percent, label=motif_ID))+
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits = c(-40, 0))+
  scale_x_discrete(position = "top")+
  theme_classic()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_blank()
  )
p_promoter

p <- plot_grid(p_promoter, p_enhancer)
p
ggsave(f_outpdf, p, width = 10, height = 4)

df = rbind(df.promoter %>% transform(Group="Promoter"), df.enhancer %>% transform(Group="Enhacner"))
write_tsv(df, f_outtxt)
