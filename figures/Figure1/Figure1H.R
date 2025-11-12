rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

plot_target_region <- function(df.id, df.bed, df.merge, chrom, region_start, region_end){
  target_region_ID0 = min(df.id$ID[df.id$Chrom==chrom])
  df.sub = df.merge %>% dplyr::filter(Chrom1==chrom, Chrom2==chrom)
  df.sub = df.sub %>% dplyr::filter(Start1>region_start & Start2>region_start) %>% dplyr::filter(End1<region_end & End2<region_end) 
  df.sub = df.sub %>% transform(diff = abs(ID1-ID2))
  
  df.bed <- df.bed %>% filter(Chrom==chrom) %>% filter(Start>region_start) %>% filter(End<region_end) %>% arrange(Start)
  df.bed$Indx <- 1:nrow(df.bed)
  df.bed = df.bed %>% transform(Start=Start/1000+target_region_ID0, End=End/1000+target_region_ID0)
  limits_min = region_start/1000 + target_region_ID0
  limits_max = region_end/1000 + target_region_ID0
  p <- ggplot()+
    geom_rect(data = df.sub, aes(xmin=ID1-0.5, xmax=ID1+0.5, ymin=ID2-0.5, ymax=ID2+0.5, fill=log2(ReadNum)))+
    geom_rect(data = df.sub, aes(xmin=ID2-0.5, xmax=ID2+0.5, ymin=ID1-0.5, ymax=ID1+0.5, fill=log2(ReadNum)))+
    geom_rect(data = df.bed, aes(xmin=Start, xmax=End, ymin=limits_min-5-Indx, ymax=limits_min+5-Indx), fill="blue")+
    geom_text(data = df.bed, aes(x=Start, y=limits_min-5-Indx, label=Indx), color="black", size=1)+
    scale_fill_gradientn(colours = c("grey80", "#FFA500", "red","darkred"))+
    geom_abline(slope = 1, intercept = 0, linetype="dashed", color="grey80")+
    labs(x="", y="") +
    theme_bw()+
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(res = list(p=p, df=df.sub))
}

f_ID = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000.ID.txt"
f_in = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000.bin.txt"
f_outtxt = "./figures/results/Figure1/Figure1H/Figure1H.At.txt"
f_outpdf = "./figures/results/Figure1/Figure1H/Figure1H.At.pdf"
chrom = "Ghir_A05" 
region_start = 11720000
region_end = 11840000

df.id = read_delim(f_ID, "\t")
df.merge = read_delim(f_in, "\t") 
f_bed = "./Supplymental/Ghir_genome/Ghir.gene.bed6"
df.bed = read_delim(f_bed, "\t", col_names = FALSE)
colnames(df.bed) <- c("Chrom", "Start", "End", "GeneName", "one", "Strand")
At_res <- plot_target_region(df.id, df.bed, df.merge, chrom, region_start, region_end)
At_res$p
ggsave(f_outpdf, At_res$p, width = 5, height = 3.5)
write_tsv(At_res$df, f_outtxt)


f_outtxt = "./figures/results/Figure1/Figure1H/Figure1H.Dt.txt"
f_outpdf = "./figures/results/Figure1/Figure1H/Figure1H.Dt.pdf"
chrom = "Ghir_D05" 
region_start = 10637000
region_end = 10738000
Dt_res <- plot_target_region(df.id, df.bed, df.merge, chrom, region_start, region_end)
Dt_res$p

ggsave(f_outpdf, Dt_res$p, width = 5, height = 3.5)
write_tsv(Dt_res$df, f_outtxt)



f_ID = "./pRIC/results/analysis/A_WT_merge/A_WT_merge_1000.ID.txt"
f_in = "./pRIC/results/analysis/A_WT_merge/A_WT_merge_1000.bin.txt"
f_outtxt = "./figures/results/Figure1/Figure1H/Figure1H.A.txt"
f_outpdf = "./figures/results/Figure1/Figure1H/Figure1H.A.pdf"

chrom = "Chr05" 
region_start = 11611000
region_end = 11712000

df.id = read_delim(f_ID, "\t")
df.merge = read_delim(f_in, "\t") 
f_bed = "./Supplymental/Gabor_genome/Gabor.gene.bed6"
df.bed = read_delim(f_bed, "\t", col_names = FALSE)
colnames(df.bed) <- c("Chrom", "Start", "End", "GeneName", "one", "Strand")

A_res <- plot_target_region(df.id, df.bed, df.merge, chrom, region_start, region_end)
A_res$p

ggsave(f_outpdf, A_res$p, width = 5, height = 3.5)
write_tsv(A_res$df, f_outtxt)
