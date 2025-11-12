rm(list = ls())
library(topGO)
library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)

load_df = function(f_protein){
  df = read_delim(f_protein, "\t")
  df$change = "nochange"
  df$change[df$padj<0.05 & df$log2FoldChange> log2(2)] <- "up"
  df$change[df$padj<0.05 & df$log2FoldChange < log2(1/2)] <- "down"
  res = list(down=df$GeneName[df$change=="down"], gene=unique(df$GeneName), changegene=df$GeneName[df$change!="nochange"])
  return(res)
}


go_data2result <- function(go_data, f_out, ontology, geneList, top_nodes=30){
  resultFisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(go_data, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(go_data, algorithm = "elim", statistic = "ks")
  allRes <- GenTable(go_data, classicFisher = resultFisher, 
                     classicKS = resultKS, elimKS = resultKS.elim, 
                     orderBy = "classicFisher", ranksOf = "classicFisher", 
                     topNodes = top_nodes)
  allRes <- allRes[allRes$Significant > 0,]
  res <- cbind(data.frame(Ontology=ontology), allRes)
  res$Gids <- sapply(res$GO.ID, function(x){return(paste(intersect(geneList, genesInTerm(go_data, x)[[1]]),collapse = " "))})
  return(res)
}

gene_GO <- function(gene_list, bg_list){
  GO_db = "./Supplymental/Ghir_genome/GO/topGo.gene.map"
  # gene_list = intersect(VIGS1_down, VIGS2_down)
  # bg_list <-  unique(df.VIGS1$GeneName)
  top_nodes = 20
  geneID2GO <- readMappings(file = GO_db)
  geneList <- factor(as.integer(bg_list %in% gene_list), levels=c(0, 1))
  names(geneList) <- bg_list
  print(sprintf("Gene list: %d, Background gene list: %d", length(gene_list), length(geneList)))
  
  MF_GOdata <- new("topGOdata", ontology="MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  BP_GOdata <- new("topGOdata", ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  CC_GOdata <- new("topGOdata", ontology="CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  MF_ref <- go_data2result(MF_GOdata, output, "MF", gene_list, 500)
  BP_ref <- go_data2result(BP_GOdata, output, "BP", gene_list, 1000)
  CC_ref <- go_data2result(CC_GOdata, output, "CC", gene_list, 500)
  df.GO <- rbind(MF_ref, BP_ref,CC_ref)
  return(df.GO)
}


f_outtxt = "./figures/results/FigureS8/FigureS8E/FigureS8E.txt"
f_outpdf = "./figures/results/FigureS8/FigureS8E/FigureS8E.pdf"

f_VIGS1_protein =  "./VIGS_RNA/results/featureCount/VIGS1.DEseq.txt"
f_VIGS2_protein =  "./VIGS_RNA/results/featureCount/VIGS2.DEseq.txt"
VIGS1_protein_gene = load_df(f_VIGS1_protein)
VIGS2_protein_gene = load_df(f_VIGS2_protein)
df.merge.GO = gene_GO(gene_list=intersect(VIGS1_protein_gene$changegene, VIGS2_protein_gene$changegene), bg_list=intersect(VIGS1_protein_gene$gene, VIGS2_protein_gene$gene))

df.BP = df.merge.GO %>% filter(Ontology=="BP")
df.BP$classicFisher = as.numeric(df.BP$classicFisher)
target_term = c("response to external stimulus", "phosphorus metabolic process", "response to nutrient levels", "response to starvation",
                "response to starvation", "cellular response to phosphate starvatio...", "response to external stimulus", "dephosphorylation")

df.BP = df.BP %>% transform(Group=ifelse(classicFisher<0.05, "anno", "none"))
df.BP$Group[df.BP$Group=="anno" & df.BP$Term %in% target_term] <- "label"
p <- ggplot(data = df.BP, aes(x = -log10(classicFisher), y=Significant, color=Group))+
  geom_point(size=0.8)+
  scale_color_manual(values = c("red",  "blue", "grey70"))+
  scale_y_continuous(limits = c(0, 100))+
  geom_text(data = df.BP[df.BP$Group=="label",], aes(x = -log10(classicFisher), y=Significant, label=Term)) +
  labs(x="-log10(p-value)", y="Gene Count")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p

ggsave(f_outpdf, p, width = 3.9, height = 3.8)
write_tsv(df.BP, f_outtxt)
