rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)

f_outtxt = "./figures/results/Figure3/Figure3G/Figure3G.txt"
f_outpdf = "./figures/results/Figure3/Figure3G/Figure3G.pdf"

f_AD = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.fisher.txt"
df.AD = read_delim(f_AD, "\t") %>% transform(Sample="Tetraploid")

class_EP <- function(x){
  chrom1 = strsplit(x[2], split = "[:]")[[1]][1]
  chrom2 = strsplit(x[4], split = "[:]")[[1]][1]
  if(chrom1==chrom2 & startsWith(chrom1, "Ghir_A")){return("intra-chrom(At)")}
  if(chrom1==chrom2 & startsWith(chrom1, "Ghir_D")){return("intra-chrom(Dt)")}
  if(chrom1!=chrom2 & startsWith(chrom1, "Ghir_D") & startsWith(chrom2, "Ghir_D")){return("inter-chrom(Dt)")}
  if(chrom1!=chrom2 & startsWith(chrom1, "Ghir_A") & startsWith(chrom2, "Ghir_A")){return("inter-chrom(At)")}
  if(chrom1!=chrom2 & startsWith(chrom1, "Ghir_A") & startsWith(chrom2, "Ghir_D")){return("At-Dt")}
  if(chrom1!=chrom2 & startsWith(chrom1, "Ghir_D") & startsWith(chrom2, "Ghir_A")){return("At-Dt")}
  return("Scaffold-")
}

df.AD$Group = apply(df.AD, 1, class_EP) 
df.AD = df.AD %>% filter(Group != "Scaffold-")
df.AD.count <- df.AD %>% group_by(Group) %>% summarise(count = n()) %>% transform(total = nrow(df.AD)) %>%
  transform(percent = count/total, Species = "Tetraploid") 

df.AD.count$Group = factor(df.AD.count$Group, levels = c("inter-chrom(At)", "intra-chrom(At)", "intra-chrom(Dt)",
                                                         "inter-chrom(Dt)", "At-Dt", "Scaffold-"))
p_AD <- ggplot() +
  geom_bar(data = df.AD.count, aes(x=Group, y=percent, fill=Group), stat = "identity")+
  geom_text(data = df.AD.count, aes(x=Group, y=percent, label=count))+
  facet_wrap(vars(Species), scales = "free")+
  labs(x="")+
  scale_fill_manual(values = c('#E27B0CFF','#FED105FF','#C190B2FF', '#6CA167FF','#145A76FF', "#1E1719FF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust=1, vjust=1))
p_AD


f_AA = "./pRIC/results/EP/A_WT_merge/A_WT_merge.EP.fisher.txt"
df.AA = read_delim(f_AA, "\t") %>% transform(Sample="Diploid")
class_AA_EP <- function(x){
  chrom1 = strsplit(x[2], split = "[:]")[[1]][1]
  chrom2 = strsplit(x[4], split = "[:]")[[1]][1]
  if(chrom1==chrom2){return("intra-chrom")}
  if(chrom1!=chrom2){return("inter-chrom")}
  return("Scaffold-")
}

df.AA$Group = apply(df.AA, 1, class_AA_EP) 
df.AA.count <- df.AA %>% group_by(Group) %>% summarise(count = n()) %>% transform(total = nrow(df.AA)) %>%
  transform(percent = count/total, Species = "Diploid") 
df.AA.count$Group = factor(df.AA.count$Group, levels = c("inter-chrom", "intra-chrom"))
p_AA <- ggplot() +
  geom_bar(data = df.AA.count, aes(x=Group, y=percent, fill=Group), stat = "identity")+
  geom_text(data = df.AA.count, aes(x=Group, y=percent, label=count))+
  facet_wrap(vars(Species), scales = "free")+
  labs(x="")+
  scale_fill_manual(values = c('#E27B0CFF','#FED105FF','#C190B2FF', '#6CA167FF','#145A76FF', "#1E1719FF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust=1, vjust=1))
p_AA

design <- "1122222"
p <- p_AA + p_AD +plot_layout(design = design)
p
df = rbind(df.AA.count, df.AD.count)
ggsave(f_outpdf, p, width = 5, height = 3)
write_tsv(df, f_outtxt)
