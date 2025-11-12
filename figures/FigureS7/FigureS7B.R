rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ggsignif)
load_chimeric <- function(f_in){
  df = read_delim(f_in, "\t") %>% filter(Group=="E-P")
  df1 = df[startsWith(df$AcceptorFeature, "pro"),]
  df1 = df1 %>% select(c("DonorFeature", "DonorChrom", "DonorCJS", "AcceptorFeature", "AcceptorChrom", "AcceptorCJS", "ReadNum"))
  colnames(df1) <- c("Enhancer","EnhancerChrom", "EnhancerSite","Promoter","PromoterChrom","PromoterSite", "ReadNum")
  df2 = df[!startsWith(df$AcceptorFeature, "pro"),]
  df2 = df2 %>% select(c("AcceptorFeature", "AcceptorChrom", "AcceptorCJS", "DonorFeature", "DonorChrom", "DonorCJS", "ReadNum"))
  colnames(df2) <- c("Enhancer","EnhancerChrom", "EnhancerSite","Promoter","PromoterChrom","PromoterSite", "ReadNum")
  dat = rbind(df1, df2)
  dat = dat %>% group_by(Enhancer,EnhancerChrom, EnhancerSite,Promoter,PromoterChrom,PromoterSite) %>% summarise(Reads = sum(ReadNum))
  return(dat)
}

f_outtxt = "./results/FigureS7/FigureS7B/FigureS7B.txt"
f_outpdf = "./results/FigureS7/FigureS7B/FigureS7B.pdf"
f_chimeric = "./pRIC/results/EP/AD_pCp_merge/AD_pCp_merge.EP.chimeric.txt"
f_MFE = "./pRIC/results/EP/AD_pCp_merge/EP_flank_MFE.txt"
df.MFE = read_delim(f_MFE, "\t") %>% transform(Sample="fl")
df.reads = load_chimeric(f_chimeric)
df = left_join(df.MFE, df.reads)

df = df %>% filter(Group=="EP")
df$ReadsGroup = ""
df$ReadsGroup[df$Reads>0 & df$Reads<=3] = "<=3"
df$ReadsGroup[df$Reads>3 & df$Reads<=5] = "4-5"
df$ReadsGroup[df$Reads>5 & df$Reads<=10] = "6-10"
df$ReadsGroup[df$Reads>10] = ">10"

df$ReadsGroup = factor(df$ReadsGroup, levels = c("<=3", "4-5", "6-10", ">10"))
p <- ggplot(data = df, aes(x=ReadsGroup, y=MFE, fill=ReadsGroup))+
  stat_boxplot(geom = "errorbar", width=0.3, position=position_dodge(0.7)) +
  # geom_violin(width=0.7, position=position_dodge(0.7))+
  geom_boxplot(outlier.shape = NA, width=0.6, position=position_dodge(0.7))+
  geom_signif(comparisons = list(c("<=3", "4-5"), c("6-10", "4-5"), c("6-10", ">10")),  
              test = t.test, y_position = c(-25, -30, -35),
              map_signif_level = F, color = "black",textsize = 4,
              tip_length=0.001, step_increase=0.004, margin_top = 0.05) +
  scale_y_continuous(limits = c(-50, 10))+
  scale_fill_brewer(palette = "Paired2")+
  labs(x="", y="MFE")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p

ggsave(f_outpdf, p, width = 3, height = 4)
write_tsv(df, f_outtxt)
