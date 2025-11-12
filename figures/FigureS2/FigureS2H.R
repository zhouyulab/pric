rm(list = ls())
library(argparse)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(reshape2)

record_miRNAfold <- function(f_miRNAfold){
  df.miRNA <- read_delim(f_miRNAfold, "\t")
  names(df.miRNA) <- c("Tid", "Chrom", "Start", "End", "Strand","PairSequence", "Pair")
  df.miRNA$loop_L <- df.miRNA$Start + sapply(strsplit(df.miRNA$Pair, ""), function(x){return(max(which(x=="(")))})
  df.miRNA$loop_R <- df.miRNA$Start + sapply(strsplit(df.miRNA$Pair, ""), function(x){return(min(which(x==")")))})
  df.miRNA$Loop_Mid <- (df.miRNA$loop_L + df.miRNA$loop_R) / 2
  df.miRNA$Length = df.miRNA$End - df.miRNA$Start 
  return(df.miRNA)
}


record_miRNA_reads <- function(miRNA_file, df.miRNA){
  df_read = read_delim(miRNA_file, "\t")[, c(1:9)]
  colnames(df_read) = c("DonorChrom", "DonorStart", "DonorStrand", "Donorblocks", "AcceptChrom", "AcceptStart", "AcceptStrand", "Acceptblocks", "ReadNum")
  df_read = df_read %>% filter(abs(DonorStart-AcceptStart)>5)
  miRNA_reads_li <- list()
  for(tid in unique(df.miRNA$Tid)){
    tid_chrom = df.miRNA$Chrom[df.miRNA$Tid==tid]
    tid_start = df.miRNA$Start[df.miRNA$Tid==tid]
    tid_end = df.miRNA$End[df.miRNA$Tid==tid]
    tid_df = df_read %>% filter(DonorChrom==tid_chrom, AcceptChrom==tid_chrom) %>% 
      filter(DonorStart>= tid_start, DonorStart<=tid_end) %>% 
      filter(AcceptStart>= tid_start, AcceptStart<=tid_end)
    if(nrow(tid_df)==0) next
    tid_df <- tid_df %>% 
      separate(col = "Donorblocks", into = c("Donorblocks",  "SDonorblocks"), remove = TRUE, sep = "[,]") %>% dplyr::select(-"SDonorblocks")%>% 
      separate(col = "Acceptblocks", into = c("Acceptblocks",  "SAcceptblocks"), remove = TRUE, sep = "[,]") %>% dplyr::select(-"SAcceptblocks")%>% 
      separate(col = "Donorblocks", into = c("DS",  "DE"), remove = FALSE, sep = "[-]") %>% 
      separate(col = "Acceptblocks", into = c("AS",  "AE"), remove = FALSE, sep = "[-]") %>%
      mutate(DS = as.numeric(DS), DE = as.numeric(DE), AS = as.numeric(AS), AE = as.numeric(AE)) 
    con1 = (tid_df$AS < tid_df$DE) & (tid_df$DS < tid_df$AS)
    con2 = (tid_df$DS < tid_df$AE) & (tid_df$AS < tid_df$DS)
    con = con1 | con2
    tid_df <- tid_df[!con,] %>% dplyr::select(-c("DS","DE","AS","AE"))
    if(nrow(tid_df)==0) next
    miRNA_reads_li[[tid]] <- tid_df
  }
  return(miRNA_reads_li)
}

block2pos <- function(block_str, start, end, strand, cnt){
  genom_pos <- do.call(c, lapply(strsplit(block_str, ";")[[1]], function(x){
    block <- as.integer(strsplit(x, "-")[[1]])
    return(block[1]:(block[2]-1))
  }))
  genom_pos <- genom_pos[genom_pos>=start & genom_pos<=(end-1)]
  if(strand=="+"){
    pos <- genom_pos-start
  }else{
    pos <- (end-1)-genom_pos
  }
  if(length(pos)==0){
    return(data.frame())
  }
  df <- data.frame(Position=pos, ReadNum=cnt)
  return(df)
}
compute_read_map_worker <- function(tid, miRNA_reads_li, miRNA_str_info){
  cat(tid, "\n")
  tmp_reads <- miRNA_reads_li[[tid]]
  tmp_reads = na.omit(tmp_reads)
  start <- miRNA_str_info$Start[miRNA_str_info$Tid==tid][1]
  end <- miRNA_str_info$End[miRNA_str_info$Tid==tid][1]
  strand <- miRNA_str_info$Strand[miRNA_str_info$Tid==tid][1]
  len <- end - start
  loop_mid <- miRNA_str_info$Loop_Mid[miRNA_str_info$Tid==tid][1]
  total_read_num <- sum(tmp_reads$ReadNum)
  
  donor_block_cnt <- do.call(rbind, lapply(1:nrow(tmp_reads), function(x){return(block2pos(tmp_reads$Donorblocks[x], start, end, strand, tmp_reads$ReadNum[x]))}))
  acceprot_block_cnt <- do.call(rbind, lapply(1:nrow(tmp_reads), function(x){return(block2pos(tmp_reads$Acceptblocks[x], start, end, strand, tmp_reads$ReadNum[x]))}))
  block_cnt_df <- rbind(donor_block_cnt, acceprot_block_cnt)
  block_cnt_df <- block_cnt_df %>% group_by(Position) %>% summarise(ReadNum=sum(ReadNum))
  block_cnt_df <- left_join(data.frame(Position=0:(end-start-1)), block_cnt_df)
  block_cnt_df$ReadNum[is.na(block_cnt_df$ReadNum)] <- 0
  block_cnt_df$Tid <- tid
  block_cnt_df$Len <- len
  block_cnt_df$PlotPos <- block_cnt_df$Position + start - loop_mid
  block_cnt_df$ReadRatio <- block_cnt_df$ReadNum / total_read_num
  block_cnt_df$Type <- "Blocks"
  
  return(block_cnt_df)
}

plot_df <- do.call(rbind, lapply(names(miRNA_reads_li), function(x){
  return(compute_read_map_worker(x, miRNA_reads_li, df.miRNA))
}))


miRNA_file = "./pRIC/results/QC/premiRNA/AD_pCp_merge.miRNA.chimeric.txt"
f_miRNAfold = "./pRIC/results/QC/premiRNA/miRNA.merge.tsv"
f_outtxt = "./figures/results/FigureS2/FigureS2H/FigureS2H.txt"
f_outpdf = "./figures/results/FigureS2/FigureS2H/FigureS2H.pdf"

df.miRNA = record_miRNAfold(f_miRNAfold)
df.miRNA.Info = df.miRNA %>% select(c("Tid", "Chrom", "Start", "End", "Strand"))
df.miRNA.Info = df.miRNA.Info[!duplicated(df.miRNA.Info),]
miRNA_reads_li <- record_miRNA_reads(miRNA_file, df.miRNA.Info)

f_anno = "./scRNA/premiRNA.ano.txt"
df.anno = read_delim(f_anno, "\t", col_names = FALSE)
f_rmDum = "./pRIC/results/QC/premiRNA/miRNA.rmDup.tsv"
df.rmDum = read_delim(f_rmDum, "\t")

miRNA_reads_li = miRNA_reads_li[names(miRNA_reads_li) %in% df.anno$X1]
miRNA_reads_li = miRNA_reads_li[names(miRNA_reads_li) %in% df.rmDum$Name]

plot_df$ReadRatio[plot_df$ReadRatio>1] <- 1
plot_df = plot_df %>% filter(Len <= 200)
df.indx = plot_df %>% select("Tid", "Len") %>% arrange(desc(Len))
df.indx = df.indx[!duplicated(df.indx),]
df.indx$Indx = 1:nrow(df.indx)
df.block <- inner_join(plot_df, df.indx)
p <- ggplot() +
  geom_rect(data=df.block, mapping = aes(xmin=PlotPos-0.5, xmax=PlotPos+0.5, ymin=Indx-0.5, ymax=Indx+0.5, fill=ReadRatio)) +
  geom_vline(xintercept = 0, size=0.2, lty=2) +
  scale_fill_gradient2(low="blue", mid="blue", high="red", midpoint = 0.3, limits=c(0, 1), breaks=c(0, 0.5, 1), labels=c("0%", "50%", "100%")) +
  labs(x="Position of the midpoint of the loop", fill="%Reads") +
  theme_bw()+
  theme(
    text = element_text(family="ArialMT", color = "black", size = 5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(family="ArialMT", color = "black", size = 5)
  )
p

ggsave(f_outpdf, p, width = 9, height = 5, units = "cm")
write_tsv(df.block, f_outtxt)


