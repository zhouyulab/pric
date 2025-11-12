rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

load_file <- function(f_in, samp){
  data <- read_delim(f_in, col_names = FALSE) 
  colnames(data) <- c("GeneID", "Strand", paste("pos",c(1:220),sep = ""))
  df.pos = data %>% filter(Strand=="template") %>% melt(measure.vars=paste("pos",c(1:220),sep = ""), variable.name = "Position", value.name = "value")  %>%
    group_by(Position) %>% dplyr::summarise(mean_sig = mean(value), sum_sig=sum(value), lower95 = Rmisc::CI(value)[3],upper95 = Rmisc::CI(value)[1]) %>%
    transform(Strand="template")
  df.neg = data %>% filter(Strand=="antisense") %>% melt(measure.vars=paste("pos",c(1:220),sep = ""), variable.name = "Position", value.name = "value")  %>%
    group_by(Position) %>% dplyr::summarise(mean_sig = mean(value), sum_sig=sum(value), lower95 = Rmisc::CI(value)[3],upper95 = Rmisc::CI(value)[1]) %>%
    transform(Strand="antisense")
  df = rbind(df.pos, df.neg) %>% transform(Sample=samp)
  return(df)
}

df.AD.rep1 = load_file("./GRO/results/meta/PRO_AD_WT_Rep1.meta.csv", "AD_WT_Rep1")
df.AD.rep2 = load_file("./GRO/results/meta/PRO_AD_WT_Rep2.meta.csv", "AD_WT_Rep2")
f_outtxt = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AD.txt"
f_outpdf = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AD.pdf"
f_outpdf_sub = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AD.sub.pdf"

df = df.AD.rep1
df$Sample = factor(df$Sample, levels = c("AD_WT_Rep2"), labels = c("Tetraploid(WT)"))
df <- df %>% separate(col = "Position", into = c("Pos", "index"), sep="[s]") %>% dplyr::select(-"Pos") %>% transform(index=as.numeric(index))
p <- ggplot(data = df, aes(x=index, y=mean_sig, color=Strand))+
  geom_line(linewidth=0.8)+
  geom_vline(xintercept = c(50), linetype="dashed", color="grey80")+
  scale_x_continuous(breaks = c(0, 50, 170, 220), labels = c("-3kb", "TSS", "TES", "3kb"))+
  labs(x="", y="")+
  geom_hline(yintercept = c(0), linetype="dashed", color="grey50")+
  scale_color_manual(values = c("#228B22","#9932CC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(color = "black"))
p
ggsave(f_outpdf, p, width = 3,height = 2.8)

p <- ggplot(data = df[df$Strand=="antisense",], aes(x=index, y=mean_sig, color=Strand))+
  geom_line(linewidth=0.8)+
  facet_wrap(vars(Sample), scales = "free_y", ncol = 3)+
  geom_vline(xintercept = c(50), linetype="dashed", color="grey80")+
  scale_x_continuous(breaks = c(0, 50, 170, 220), labels = c("-3kb", "TSS", "TES", "3kb"))+
  labs(x="", y="")+
  geom_hline(yintercept = c(0), linetype="dashed", color="grey50")+
  scale_color_manual(values = c("#228B22","#9932CC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(color = "black"))
p
ggsave(f_outpdf_sub, p, width = 3,height = 2.8)
write_tsv(df, f_outtxt)


df.AA.rep1 = load_file("./GRO/results/meta/PRO_AA_WT_Rep1.meta.csv", "AA_WT_Rep1")
df.AA.rep2 = load_file("./GRO/results/meta/PRO_AA_WT_Rep2.meta.csv", "AA_WT_Rep2")
f_outtxt = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AA.txt"
f_outpdf = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AA.pdf"
f_outpdf_sub = "./github/figures/results/FigureS3/FigureS3E/Figure3E.AA.sub.pdf"

df = df.AA.rep1
df$Sample = factor(df$Sample, levels = c("AD_WT_Rep2"), labels = c("Diploid(WT)"))
df <- df %>% separate(col = "Position", into = c("Pos", "index"), sep="[s]") %>% dplyr::select(-"Pos") %>% transform(index=as.numeric(index))
p <- ggplot(data = df, aes(x=index, y=mean_sig, color=Strand))+
  geom_line(linewidth=0.8)+
  geom_vline(xintercept = c(50), linetype="dashed", color="grey80")+
  scale_x_continuous(breaks = c(0, 50, 170, 220), labels = c("-3kb", "TSS", "TES", "3kb"))+
  labs(x="", y="")+
  geom_hline(yintercept = c(0), linetype="dashed", color="grey50")+
  scale_color_manual(values = c("#228B22","#9932CC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(color = "black"))
p
ggsave(f_outpdf, p, width = 3,height = 2.8)

p <- ggplot(data = df[df$Strand=="antisense",], aes(x=index, y=mean_sig, color=Strand))+
  geom_line(linewidth=0.8)+
  facet_wrap(vars(Sample), scales = "free_y", ncol = 3)+
  geom_vline(xintercept = c(50), linetype="dashed", color="grey80")+
  scale_x_continuous(breaks = c(0, 50, 170, 220), labels = c("-3kb", "TSS", "TES", "3kb"))+
  labs(x="", y="")+
  geom_hline(yintercept = c(0), linetype="dashed", color="grey50")+
  scale_color_manual(values = c("#228B22","#9932CC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(color = "black"))
p
ggsave(f_outpdf_sub, p, width = 3,height = 2.8)
write_tsv(df, f_outtxt)
