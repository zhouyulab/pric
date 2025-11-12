rm(list = ls())
library(readr)
library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(corrplot)

f_outtxt = "./figures/results/FigureS2/FigureS2M/FigureS2M.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2M/FigureS2M.pdf"

f_2023_NG = "./GWAS/2023_NG_eQTL/download/Regulatory_controls_of_duplicated_gene.xlsx"
df.At2Dt.map = read_excel(f_2023_NG, sheet = "Supp Table 4", skip = 1)
df.At2Dt.map = df.At2Dt.map %>% select(c("Gene in the At subgenome", "Gene in the Dt subgenome"))
colnames(df.At2Dt.map) <- c("AtGene", "DtGene")

f_AD = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge.intra.Chimeric.txt"
df.AD = read_delim(f_AD, "\t")
df.AD = df.AD %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = '[|]')
AD_total_reads = sum(df.AD$ReadNum)
df.AD = df.AD %>% transform(NormRead = ReadNum*1000000/AD_total_reads)
df.AD = df.AD %>% select(c("GeneName","ReadNum","NormRead"))

df.AD.A = df.AD %>% filter(GeneName %in% df.At2Dt.map$AtGene)
colnames(df.AD.A) <- c("AtGene", "AtReadNum", "AtNormRead")
df.AD.D = df.AD %>% filter(GeneName %in% df.At2Dt.map$DtGene)
colnames(df.AD.D) <- c("DtGene", "DtReadNum", "DtNormRead")

f_A = "./pRIC/results/analysis/A_WT_merge.intra.Chimeric.txt"
df.A = read_delim(f_A, "\t")
df.A = df.A %>% separate(col = "GeneName", into = c("GeneType", "GeneName"), sep = '[|]')
A_total_reads = sum(df.A$ReadNum)
df.A = df.A %>% transform(NormRead = ReadNum*1000000/A_total_reads)
df.A = df.A %>% filter(GeneType=="protein") %>% dplyr::select(c("GeneName","ReadNum","NormRead"))

f_gene_map = "./Supplymental/Gabor_genome/anno/annotation/At/gene_map.txt"
df.gene.map = read_delim(f_gene_map, col_names = FALSE)

colnames(df.gene.map) <- c("AGene", "AtGene", "isFirst")
colnames(df.A) <- c("AGene", "AReadNum", "ANormRead")
df.A = inner_join(df.A, df.gene.map)

df = inner_join(df.AD.A, df.At2Dt.map)
df = inner_join(df, df.AD.D)
df = df %>% transform(AGene=NA, AReadNum=NA, ANormRead=NA)
df.add = data.frame()
for(i in 1:nrow(df)){
  At_gene = df$AtGene[i]
  Dt_gene = df$DtGene[i]
  if(!(At_gene %in% df.A$AtGene) & !(Dt_gene %in% df.A$AtGene)) next()
  if((At_gene %in% df.A$AtGene) & !(Dt_gene %in% df.A$AtGene)){
    df$AGene[i] = df.A$AGene[df.A$AtGene==At_gene]
    df$AReadNum[i] = df.A$AReadNum[df.A$AtGene==At_gene]
    df$ANormRead[i] = df.A$ANormRead[df.A$AtGene==At_gene]
  }
  if(!(At_gene %in% df.A$AtGene) & (Dt_gene %in% df.A$AtGene)){
    df$AGene[i] = df.A$AGene[df.A$AtGene==Dt_gene]
    df$AReadNum[i] = df.A$AReadNum[df.A$AtGene==Dt_gene]
    df$ANormRead[i] = df.A$ANormRead[df.A$AtGene==Dt_gene]
  }
  if((At_gene %in% df.A$AtGene) & (Dt_gene %in% df.A$AtGene)){
    aGene_At = df.A$AGene[df.A$AtGene==At_gene]
    aGene_Dt = df.A$AGene[df.A$AtGene==Dt_gene]
    common_gene = intersect(aGene_At, aGene_Dt)
    if(length(common_gene)==1){
      df$AGene[i] = df.A$AGene[df.A$AtGene==common_gene][1]
      df$AReadNum[i] = df.A$AReadNum[df.A$AtGene==common_gene][1]
      df$ANormRead[i] = df.A$ANormRead[df.A$AtGene==common_gene][1]
    }else{
      a_gene = df.A$AGene[df.A$AtGene %in% c(At_gene, Dt_gene)]
      a_readNum = df.A$AReadNum[df.A$AtGene %in% c(At_gene, Dt_gene)]
      a_normRead = df.A$ANormRead[df.A$AtGene %in% c(At_gene, Dt_gene)]
      temp_df = data.frame(AtGene=df$AtGene[i], AtReadNum=df$AtReadNum[i], AtNormRead=df$AtNormRead[i], DtGene=df$DtGene[i], DtReadNum=df$DtReadNum[i],
                           DtNormRead=df$DtNormRead[i], AGene=a_gene, AReadNum=a_readNum, ANormRead=a_normRead)
      df.add = rbind(df.add, temp_df)
    }
  }
  
}

df.filter = rbind(na.omit(df), df.add)
dat = df.filter %>% select(c("ANormRead", "AtNormRead", "DtNormRead"))
dat = cor(log2(dat))
COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)

pdf(f_outpdf, width = 5, height = 5)
corrplot(corr=dat, col = COL2('PiYG'), type="upper",method="ellipse",diag = FALSE,
         tl.col="black", addCoef.col = "black",rect.col = 'black')
dev.off()

write_tsv(as.data.frame(dat), f_outtxt)






