rm(list = ls())
library(grid)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

plot_GRO2 <- function(f_in){
  df = read_delim(f_in, "\t")
  
  df = df %>% transform(Anno = ifelse(Position=="Intergenic", "Intergenic", "Genebody"))
  df.pos = df %>% group_by(Anno) %>% summarise(num=n())
  df.pos$Anno = factor(df.pos$Anno, levels = c("Intergenic", "Genebody"))
  df.pos = df.pos %>% arrange(desc(Anno))
  df.pos = df.pos  %>% transform(Total=nrow(df)) %>% transform(percent=num/Total) %>%
    transform(label=paste(Anno, "(", num, ";", round(percent*100, 2), "%", ")", sep = ""))
  df.pos$End = cumsum(df.pos$num)
  df.pos$Start = c(1, df.pos$End[1:(nrow(df.pos)-1)])
  
  df.gene = df %>% filter(Anno=="Genebody") 
  df.gene$Position = factor(df.gene$Position, levels = c("upstream", "gene","downstream"))
  df.gene.pos = df.gene %>% transform(Total=nrow(df.gene)) %>% group_by(Position, Total) %>% summarise(num=n()) %>% arrange(Position) %>% transform(End=cumsum(num))
  df.gene.pos = df.gene.pos %>% transform(percent=num/Total) %>%
    transform(label=paste(Position, "(", num, ";", round(percent*100, 2), "%", ")", sep = ""))
  df.gene.pos$Start = c(1, df.gene.pos$End[1:(nrow(df.gene.pos)-1)])
  
  p <- ggplot()+
    geom_rect(data = df.gene.pos, aes(xmin = 0.25, xmax = 1.5+0.3,ymin = Start, ymax = End, fill = label), linewidth = 0.2)+
    geom_rect(data = df.pos, aes(xmin = 0.5, xmax = 1.3, ymin = Start, ymax = End, fill = label), alpha = 1, color = "white", linewidth = 0.42)+
    annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "white")+
    scale_size_identity() +
    scale_color_identity() +
    scale_fill_manual(values = c("#CCCC99","#FFCCFF", "#CCCCFF","#FF9933FF","#62AFD7FF","#9999FFFF", "#66CC33FF"))+
    coord_polar(theta = "y", clip = "off", start = 10) +
    labs(title =NULL) +
    theme_void() +
    theme(plot.background = element_rect(color = "white", fill = "white"),
          text = element_text(color ="black"),
          legend.title = element_blank(),
          plot.margin = margin(t = 0, b = 0, l = 4, r = 4))
  p
  
  outdf = rbind(df.gene.pos %>% dplyr::select(c("Position", "num", "Total", "percent", "label")) %>% dplyr::rename("Anno"="Position"),
                df.pos %>% dplyr::select(c("Anno", "num", "Total", "percent", "label")))
  return(res = list(p=p, df=outdf))
}

f_AD = "./GRO/results/feature/AD_GRO.filter.anno.bed"
f_AA = "./GRO/results/feature/AA_GRO.filter.anno.bed"
f_outtxt = "./figures/results/Figure2/Figure2A/Figure2A.txt"
f_outpdf = "./figures/results/Figure2/Figure2A/Figure2A.pdf"

AD_res = plot_GRO2(f_AD)
AA_res = plot_GRO2(f_AA)
p <- plot_grid(AA_res$p, AD_res$p, labels = c("AA", "AD"))
p

df = rbind(AD_res$df %>% transform(Species="Tetraploid"),
           AA_res$df %>% transform(Species="Diploid"))
ggsave(f_outpdf, p, width = 8, height = 3)
write_tsv(df, f_outtxt)



