rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

f_ID = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000.ID.txt"
f_pCp = "./pRIC/results/analysis/AD_pCp_merge/AD_pCp_merge_1000.bin.txt"
f_outtxt = "./figures/results/Figure1/Figure1F/Figure1F.txt"
f_outpdf = "./figures/results/Figure1/Figure1F/Figure1F.pdf"

df.id = read_delim(f_ID, "\t")
df.merge = read_delim(f_pCp, "\t")
f_bed = "./Supplymental/Ghir_genome/Ghir.gene.bed6"
df.bed0 = read_delim(f_bed, "\t", col_names = FALSE)
colnames(df.bed0) <- c("Chrom", "Start", "End", "GeneName", "one", "Strand")

chrom = "Ghir_D04" 
region_start = 11060*1000
region_end = 11110*1000

target_region_ID0 = min(df.id$ID[df.id$Chrom==chrom])
df.sub = df.merge %>% filter(Chrom1==chrom, Chrom2==chrom)
df.sub = df.sub %>% filter(Start1>region_start & Start2>region_start) %>% filter(End1<region_end & End2<region_end) 
df.sub = df.sub %>% transform(diff = abs(ID1-ID2))

df.bed <- df.bed0 %>% filter(Chrom==chrom) %>% filter(Start>region_start) %>% filter(End<region_end) %>% arrange(Start)
df.bed$Indx <- 1:nrow(df.bed)
df.bed = df.bed %>% transform(Start=Start/1000+target_region_ID0, End=End/1000+target_region_ID0)
limits_min = region_start/1000 + target_region_ID0
limits_max = region_end/1000 + target_region_ID0
df.sub$ReadNum[df.sub$ReadNum>100] = 100
p <- ggplot()+
  geom_rect(data = df.sub, aes(xmin=ID1-0.5, xmax=ID1+0.5, ymin=ID2-0.5, ymax=ID2+0.5, fill=log2(ReadNum)))+
  geom_rect(data = df.sub, aes(xmin=ID2-0.5, xmax=ID2+0.5, ymin=ID1-0.5, ymax=ID1+0.5, fill=log2(ReadNum)))+
  geom_rect(data = df.bed, aes(xmin=Start, xmax=End, ymin=limits_min-5-Indx, ymax=limits_min+5-Indx), fill="blue")+
  geom_text(data = df.bed, aes(x=Start, y=limits_min-5-Indx, label=Indx), color="black", size=1)+
  scale_fill_gradientn(colours = c("grey80", "#FFA500", "red","darkred"))+
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color="grey80")+
  scale_x_continuous( limits = c(limits_min, limits_max))+
  scale_y_continuous( limits = c(limits_min-10, limits_max))+
  labs(x="", y="") +
  theme_bw()+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p

ggsave(f_outpdf, p, width = 5, height = 3.5)
write_tsv(df.sub, f_outtxt)
