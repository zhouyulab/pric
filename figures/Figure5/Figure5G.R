rm(list = ls())
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(topGO)
library(ggsignif)

load_EP_gene <- function(f_in, cutoff=0){
  EP_gene = c()
  df = read_delim(f_in, "\t") %>% filter(Group=="E-P") %>% filter(Reads>cutoff)
  for(i in 1:nrow(df)){
    gene = gsub("pro", "", df$Promoter[i])
    EP_gene = c(EP_gene, gene)
  }
  return(EP_gene)
}

f_outtxt = "./figures/results/Figure5/Figure5G/Figure5G.txt"
f_outpdf = "./figures/results/Figure5/Figure5G/Figure5G.pdf"

f_WT = "./pRIC/results/EP/AD_WT_merge/AD_WT_merge.EP.fisher.txt"
f_fl = "./pRIC/results/EP/AD_fl_merge/AD_fl_merge.EP.fisher.txt"
f_RNA = "./AD_RNA/results/featureCount/RNA.FPKM.txt"

df_RNA = read_delim(f_RNA, "\t")
df_RNA = df_RNA %>% transform(WT = (WT_Rep1+WT_Rep2)/2, fl = (fl_Rep1+fl_Rep2)/2)

WT_gene = load_EP_gene(f_WT, cutoff = 3)
fl_gene = load_EP_gene(f_fl, cutoff = 3)
common_gene = intersect(WT_gene, fl_gene)

f_DEseq = "./AD_RNA/results/featureCount/RNA.DEseq.res"
df.DEseq = read_delim(f_DEseq, "\t")
df.DEseq$Group = "None"
df.DEseq$Group[df.DEseq$GeneName %in% common_gene] <- "common"
df.DEseq$Group[df.DEseq$Group=="None" & df.DEseq$GeneName %in% WT_gene] <- "WT-specific"
df.DEseq$Group[df.DEseq$Group=="None" & df.DEseq$GeneName %in% fl_gene] <- "fl-specific"
df.DEseq = df.DEseq %>% filter(Group != "None")
df.DEseq$Group = factor(df.DEseq$Group, levels = c("WT-specific", "common", "fl-specific"))

GO_db = "./Supplymental/Ghir_genome/GO/topGo.gene.map"
gene_list <- unique(df.DEseq$GeneName[df.DEseq$Group=="WT-specific"])
bg_list <-  unique(df.DEseq$GeneName[df.DEseq$Group %in% c("WT-specific","fl-specific")])

top_nodes = 10
geneID2GO <- readMappings(file = GO_db)
geneList <- factor(as.integer(bg_list %in% gene_list), levels=c(0, 1))
names(geneList) <- bg_list
print(sprintf("Gene list: %d, Background gene list: %d", length(gene_list), length(geneList)))
BP_GOdata <- new("topGOdata", ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
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
BP_ref <- go_data2result(BP_GOdata, output, "BP", gene_list, 1000)
target_BP <- c("plant epidermis development")
target_BP <- c("seed trichome differentiation","seed trichome elongation","trichome differentiation",
               "seed development", "plant epidermis development", "fruit development")
df.BP = BP_ref
df.BP$classicFisher = as.numeric(df.BP$classicFisher)
df.BP = df.BP %>% transform(Group=ifelse(classicFisher<0.05, "**", "NS"))
df.BP$anno = ""
df.BP$anno[df.BP$Term %in% target_BP] = df.BP$Term[df.BP$Term %in% target_BP]
df.BP$Group[df.BP$Term %in% target_BP] <- "***"
p <- ggplot(data = df.BP, aes(x=-log10(classicFisher), y=Significant, color=Group))+
  geom_point(size=0.8)+
  labs(x="-log10(p-value)", y="gene count")+
  scale_color_manual(values = c("red", "blue","grey70"))+
  geom_text(aes(x=-log10(classicFisher), y=Significant, label=anno))+
  scale_y_continuous(limits = c(0, 200))+
  geom_vline(xintercept = c(-log10(0.1)), linetype="dashed")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none")
p
ggsave(f_outpdf, p, width = 3, height = 2.8)
write_tsv(df.BP, f_outtxt)
