rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggbreak)

f_outtxt = "./figures/results/FigureS10/FigureS10F/FigureS10F.promoter.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10F/FigureS10F.promoter.pdf"
f_indir = "./plantTFDB/results/TF_enrich/"
promoter_list = list.files(f_indir, pattern = "promoter.txt")
df_li = list()
for (promoter in promoter_list){
  TF = strsplit(promoter, split = "\\.")[[1]][1] 
  f_in = file.path(f_indir, promoter)
  df_temp = read_delim(f_in) %>% select(c("Idex", "Percent")) %>% transform(TF = TF)
  df_li[[TF]] = df_temp
}

df = Reduce(rbind, df_li)
target_TF = c("C2H2", "ERF", "LBD", "GATA", "MYB", "BBR-BPC")
df = df %>% transform(label = ifelse(TF %in% target_TF, TF, "other")) %>% arrange(label)
df$label = factor(df$label, levels = c("ERF", "C2H2", "LBD", "GATA", "MYB", "BBR-BPC", "other"))
df$Group="Promoter"
p <- ggplot(data = df, aes(x=Idex, y=Percent, group=TF, color=label))+
  geom_line(linewidth=0.8)+
  facet_wrap(vars(Group))+
  labs(x="", y="Precent (%)")+
  scale_color_manual(values = c("#C70E7B", "#FC6882", "#EF7C12", "#007BC3", "#54BCD1", "#8FDA04", "grey80"))+
  scale_x_continuous(limits = c(1500, 3500),breaks = c(1500, 2500, 3500), labels = c("-1kb", "CJS", "1kb"))+
  scale_y_break(breaks = c(1.5,2),space = 0.1,scales = 0.5,expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"))
p

ggsave(f_outpdf, p, width = 4, height = 3)
outdf = df %>% dplyr::filter(Idex>= 1500, Idex <= 3500)
write_tsv(outdf, f_outtxt)


rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

f_outtxt = "./figures/results/FigureS10/FigureS10F/FigureS10F.enhancer.txt"
f_outpdf = "./figures/results/FigureS10/FigureS10F/FigureS10F.enhancer.pdf"
f_indir =  "./plantTFDB/results/TF_enrich/"
enhancer_list = list.files(f_indir, pattern = "enhancer.txt")
df_li = list()
for (enhancer in enhancer_list){
  TF = strsplit(enhancer, split = "\\.")[[1]][1] 
  f_in = file.path(f_indir, enhancer)
  df_temp = read_delim(f_in) %>% select(c("Idex", "Percent")) %>% transform(TF = TF)
  df_li[[TF]] = df_temp
}

df = Reduce(rbind, df_li)
target_TF = c("Dof", "ERF", "BBR-BPC")
df = df %>% transform(label = ifelse(TF %in% target_TF, TF, "other")) %>% arrange(label)
df$label = factor(df$label, levels = c("ERF", "C2H2", "LBD", "GATA", "MYB", "BBR-BPC", "Dof","other"))
df$Group="Enhancer"
p <- ggplot(data = df, aes(x=Idex, y=Percent, group=TF, color=label))+
  geom_line(linewidth=0.8)+
  facet_wrap(vars(Group), scales = "free")+
  labs(x="", y="Precent (%)")+
  scale_color_manual(values = c("#C70E7B", "#007BC3", "#54BCD1", "grey80"))+
  scale_x_continuous(limits = c(1500, 3500),breaks = c(1500, 2500, 3500), labels = c("-1kb", "CJS", "1kb"))+
  scale_y_continuous(limits = c(0, 1.1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"))
p

ggsave(f_outpdf, p, width = 4, height = 3)
outdf = df %>% dplyr::filter(Idex>= 1500, Idex <= 3500)
write_tsv(outdf, f_outtxt)
