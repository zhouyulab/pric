rm(list = ls())
library(topGO)
library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)

f_outtxt = "./figures/results/FigureS8/FigureS8D/FigureS8D.txt"
f_outpdf = "./figures/results/FigureS8/FigureS8D/FigureS8D.pdf"

load_df = function(f_protein){
  df = read_delim(f_protein, "\t")
  df$change = "nochange"
  df$change[df$padj<0.05 & df$log2FoldChange> log2(2)] <- "up"
  df$change[df$padj<0.05 & df$log2FoldChange < log2(1/2)] <- "down"
  res = list(down=df$GeneName[df$change=="down"], gene=unique(df$GeneName), changegene=df$GeneName[df$change!="nochange"])
  return(res)
}

plot_overlap = function(gene1_li, gene2_li, label1, label2){
  p_down <-  venn.diagram(
    list(gene1_li, gene2_li),
    category.names = c(label1 , label2),
    filename = NULL, cex = 2, cat.cex=1.5, cat.default.pos='text',
    fill = c("#FC6882", "#007BC3FF"), cat.pos = c(1,1), col = "transparent", lwd = 2,
    main = paste("VIGS down gene")
  )
  return(p_down)
}

f_VIGS1_protein =  "./VIGS_RNA/results/featureCount/VIGS1.DEseq.txt"
f_VIGS2_protein =  "./VIGS_RNA/results/featureCount/VIGS2.DEseq.txt"
f_VIGS1_lncRNA =  "./VIGS_RNA/results/lncRNA/VIGS1.DEseq.txt"
f_VIGS2_lncRNA =  "./VIGS_RNA/results/lncRNA/VIGS2.DEseq.txt"

VIGS1_protein_gene = load_df(f_VIGS1_protein)
VIGS2_protein_gene = load_df(f_VIGS2_protein)
VIGS1_lncRNA_gene = load_df(f_VIGS1_lncRNA)
VIGS2_lncRNA_gene = load_df(f_VIGS2_lncRNA)

p_change_overlap <- plot_overlap(VIGS1_protein_gene$changegene,
                                 VIGS2_protein_gene$changegene,
                                 "VIGS1", "VIGS2")
pdf(f_outpdf)
grid.draw(p_change_overlap)
dev.off()

outdf = rbind(data.frame(GeneName=VIGS1_protein_gene$changegene, Sample="VIGS1"),
              data.frame(GeneName=VIGS2_protein_gene$changegene, Sample="VIGS2"))
write_tsv(outdf, f_outtxt)
