rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

science_point <- function(pvalue, n){
  temp <- as.character(pvalue)
  if("e" %in% base::strsplit(temp, split = "")[[1]]){
    temp <- base::strsplit(temp, split = "e")[[1]]
    num <- round(as.numeric(temp[1]), n)
    index <- as.numeric(temp[2])
    # p <- bquote(.(num)^".(-index)")
    p <- num * 10^index
  }else{
    p <- round(pvalue, n)
  }
  return(p)
}

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
calculate_pvalue <- function(df){
  df.pvalue = data.frame()
  for(te in unique(df$label)){
    df.enhancer = df %>% filter(label==te, Feature=="enhancer")  %>% group_by(Total) %>% summarise(Num=round(mean(Num),0))
    df.enhancer = df.enhancer %>% transform(noNum=Total-Num)
    enhancer.mat = as.matrix(df.enhancer[, c("Num", "noNum")])
    enhancer_pvalue = chisq.test(enhancer.mat)$p.value
    
    df.promoter = df %>% filter(label==te, Feature=="promoter")  %>% group_by(Total) %>% summarise(Num=round(mean(Num),0))
    df.promoter = df.promoter %>% transform(noNum=Total-Num)
    promoter.mat = as.matrix(df.promoter[, c("Num", "noNum")])
    promoter_pvalue = chisq.test(promoter.mat)$p.value
    
    temp_df = data.frame(Feature=c("enhancer", "promoter"), pvalue=c(annotate_pvalue(enhancer_pvalue), annotate_pvalue(promoter_pvalue)), label=te)
    df.pvalue = rbind(df.pvalue, temp_df)
  }
  return(df.pvalue)
}

f_outtxt = "./figures/results/Figure4/Figure4D/Figure4D.txt"
f_outpdf = "./figures/results/Figure4/Figure4D/Figure4D.pdf"

f_indir = "./pRIC/results/TE/chimeric/"
filenames = list.files(f_indir, pattern = "*250.TE.statis.txt")
df_li = list()
for(filename in filenames){
  f_in = file.path(f_indir, filename)
  df = read_delim(f_in, "\t") %>% transform(filename=filename)
  df_li[[filename]] = df
}

df = Reduce(rbind, df_li)
df = df %>% transform(Group=ifelse(startsWith(filename, "EP"), "EP", "Random"))

RNA_TE = c("Unknown", "LTR/Gypsy","Simple_repeat", "Low_complexity", "LINE/L1", "LTR/Copia", "SINE/tRNA", "tRNA", "rRNA", "LTR/Caulimovirus", "LTR/ERV4", "Satellite")
DNA_TE = c("DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/hAT", "DNA/hAT-Tag1", "DNA/hAT-Ac", "DNA", "DNA/hAT-Tip100", "RC/Helitron", "DNA/CMC-EnSpm", "DNA/TcMar", "DNA/Maverick")
target_TE = c("Unknown", "LTR/Gypsy", "LTR/Copia", "DNA/MULE-MuDR", "DNA/hAT-Ac", "LINE/L1", "Simple_repeat", "Low_complexity")
df = df %>% transform(TE = ifelse(Transposable %in% target_TE, Transposable, "other")) %>% transform(RAN_Group=ifelse(Transposable %in% RNA_TE, "RNA", "DNA"))
df = df %>% transform(label=paste(RAN_Group, TE, sep = "-"))

df.count = df %>% group_by(label, Total, filename, Feature) %>% summarise(Num=sum(Num)) %>% transform(Percent=Num*100/Total)
df.count = df.count %>% transform(Group=ifelse(startsWith(filename, "EP"), "EP", "Random"))
df.count$label = factor(df.count$label, levels = c("RNA-Simple_repeat", "RNA-LTR/Gypsy", "RNA-LINE/L1", "RNA-LTR/Copia", "RNA-other",
                                                   "DNA-DNA/MULE-MuDR", "DNA-DNA/hAT-Ac", "DNA-other", "RNA-Unknown","RNA-Low_complexity"))
df.count$Group = factor(df.count$Group, levels = c("Random", "EP"))
df.pvalue = calculate_pvalue(df.count) %>% transform(Group="EP")
df.mean = df.count %>% group_by(Group, label, Feature) %>% summarise(meanPercent = mean(Percent)) 
p <- ggplot(data = df.count, aes(x=label, y=Percent, color=Group))+
  geom_bar(data = df.mean, aes(x=label, y=meanPercent, fill=Group), stat = "identity", position = position_dodge(), color=NA)+
  geom_text(data = df.pvalue, aes(x=label, y=16, label=pvalue))+
  facet_wrap(vars(Feature))+
  labs(x="", y="Percent (%)") +
  geom_point(position = position_dodge(width = 0.9))+
  scale_fill_brewer(palette = "Paired")+
  scale_color_manual(values = c("grey60", "black"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5))
p

ggsave(f_outpdf, p, width = 6.5, height = 3.5)
write_tsv(df.count, f_outtxt)
