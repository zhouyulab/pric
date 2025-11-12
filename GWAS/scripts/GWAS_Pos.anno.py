import argparse
import os
import sys
from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import numpy as np

class GeneLine(object):
    def __init__(self, chrom, start, end, transcriptID, gene_name, strand, pos):
        # print(line)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.transcriptID = transcriptID
        self.gene_name = gene_name
        self.strand = strand
        self.position = pos


class transcriptFeature(object):
    def __init__(self, bed_line):
       chrom, start, end, transcriptID, zero, strand, cds_start, cds_end, source, exonNum, length_li, start_li = bed_line
       self.transcriptID = transcriptID
       self.Chrom = chrom
       self.Start = int(start)
       self.End = int(end)
       self.Strand = strand
       self.cdsStart = int(cds_start)
       self.cdsEnd = int(cds_end)
       self.exonNum = int(exonNum)
       self.length_li = list(map(int, length_li.split(",")[:-1]))
       self.start_li = list(map(int, start_li.split(",")[:-1]))
       self.Region = None
       self.Blocks = None

    def get_region(self):
        self.Region = IntervalTree()
        if self.Blocks is None:
            self.get_blocks()
        for blcok_start, block_end in self.Blocks:
            assert blcok_start <= block_end
            if block_end <= self.cdsStart:
                block_anno = "UTR5" if self.Strand == "+" else "UTR3"
                self.Region.insert(int(blcok_start), int(block_end), block_anno)
            elif blcok_start <= self.cdsStart and block_end >= self.cdsStart and block_end <= self.cdsEnd:
                block_anno = "UTR5" if self.Strand == "+" else "UTR3"
                self.Region.insert(int(blcok_start), int(self.cdsStart), block_anno)
                self.Region.insert(int(self.cdsStart), int(block_end), "exon")
            elif blcok_start >= self.cdsStart and block_end <= self.cdsEnd:
                self.Region.insert(int(blcok_start), int(block_end), "exon")
            elif blcok_start >= self.cdsStart and blcok_start <= self.cdsEnd and block_end >= self.cdsEnd:
                block_anno = "UTR3" if self.Strand == "+" else "UTR5"
                self.Region.insert(int(blcok_start), int(self.cdsEnd), "exon")
                self.Region.insert(int(self.cdsEnd), int(block_end), block_anno)
            elif blcok_start >= self.cdsEnd:
                block_anno = "UTR3" if self.Strand == "+" else "UTR5"
                self.Region.insert(int(blcok_start), int(block_end), block_anno)
            elif blcok_start <= self.cdsStart and block_end >= self.cdsEnd:
                block_anno1 = "UTR5" if self.Strand == "+" else "UTR3"
                block_anno2 = "UTR3" if self.Strand == "+" else "UTR5"
                self.Region.insert(int(blcok_start), self.cdsStart, block_anno1)
                self.Region.insert(self.cdsStart, self.cdsEnd, "exon")
                self.Region.insert(self.cdsEnd, int(block_end), block_anno2)
            else:
                print(blcok_start, block_end, "anno:", self.cdsStart, self.cdsEnd)

        for i in range(self.exonNum-1):
            intron_start, intron_end = self.Blocks[i][1], self.Blocks[i+1][0]
            if intron_end <= self.cdsStart:
                block_anno = "UTR5_intron" if self.Strand == "+" else "UTR3_intron"
                self.Region.insert(intron_start, intron_end, block_anno)
            elif intron_start <= self.cdsStart and intron_end >= self.cdsStart and intron_end <= self.cdsEnd:
                block_anno = "UTR5_intron" if self.Strand == "+" else "UTR3_intron"
                self.Region.insert(intron_start, self.cdsStart, block_anno)
                self.Region.insert(self.cdsStart, intron_end, "intron")
            elif intron_start >= self.cdsStart and intron_end <= self.cdsEnd:
                self.Region.insert(intron_start, intron_end, "intron")
            elif intron_start >= self.cdsStart and intron_start <= self.cdsEnd and intron_end >= self.cdsEnd:
                block_anno = "UTR3_intron" if self.Strand == "+" else "UTR5_intron"
                self.Region.insert(intron_start, self.cdsEnd, "intron")
                self.Region.insert(self.cdsEnd, intron_end, block_anno)
            elif intron_start >= self.cdsEnd:
                block_anno = "UTR3_intron" if self.Strand == "+" else "UTR5_intron"
                self.Region.insert(intron_start, intron_end, block_anno)
            elif intron_start <= self.cdsStart and intron_end >= self.cdsEnd:
                block_anno1 = "UTR5_intron" if self.Strand == "+" else "UTR3_intron"
                block_anno2 = "UTR3_intron" if self.Strand == "+" else "UTR5_intron"
                self.Region.insert(intron_start, self.cdsStart, block_anno1)
                self.Region.insert(self.cdsStart, self.cdsEnd, "intron")
                self.Region.insert(self.cdsEnd, intron_end, block_anno2)
            else:
                print(intron_start, intron_end, "anno:", self.cdsStart, self.cdsEnd)

    def get_blocks(self):
        block_li = list()
        for i in range(self.exonNum):
            block_start = self.Start + self.start_li[i]
            block_end = block_start + self.length_li[i]
            block_li.append((block_start, block_end))
        block_li = sorted(block_li, key=lambda x: x[0], reverse=False)
        self.Blocks = block_li
        return block_li

    def get_rel_pos(self, site):
        site = int(site)
        if self.Strand  == "+":
            rel_pos = (site - self.Start) / (self.End - self.Start)
        else:
            rel_pos = (self.End - site) / (self.End - self.Start)
        return round(rel_pos, 6)


def load_transcript_info(f_bed):
    TransInfo = dict()
    with open(f_bed) as f:
        for line in f:
            bed_line = line.strip().split("\t")
            transcript_info = transcriptFeature(bed_line)
            transcript_info.get_region()
            transcriptID = bed_line[3]
            TransInfo[transcriptID] = transcript_info
    return TransInfo


def load_gene(f_bed, flank=1000):
    GeneInfo = dict()
    with open(f_bed) as f:
        for line in f:
            data = line.strip().split("\t")
            chrom, start, end, transcriptID, gene_name, strand = data[0], data[1], data[2], data[3], data[4], data[5]
            gene_line = GeneLine(chrom, start, end, transcriptID, gene_name, strand, "gene")
            start, end = int(start), int(end)
            if strand == "+":
                up_start, up_end = max(1, start-flank), start
                up_line = GeneLine(chrom, up_start, up_end, transcriptID, gene_name, strand, "upstream")
                down_start, down_end = end, end + flank
                down_line = GeneLine(chrom, down_start, down_end, transcriptID, gene_name, strand, "downstream")
            else:
                up_start, up_end = end, end + flank
                up_line = GeneLine(chrom, up_start, up_end, transcriptID, gene_name, strand, "upstream")
                down_start, down_end = max(1, start-flank), start
                down_line = GeneLine(chrom, down_start, down_end, transcriptID, gene_name, strand, "downstream")        

            key = chrom
            if key not in GeneInfo:
                GeneInfo[key] = IntervalTree()
            GeneInfo[key].insert(int(start), int(end), gene_line)
            GeneInfo[key].insert(int(up_start), int(up_end), up_line)
            GeneInfo[key].insert(int(down_start), int(down_end), down_line)
    return GeneInfo

def find_gid(chrom, start, end, GeneInfo):
    sidx = int(start)
    eidx = int(end)
    key = chrom
    if key not in GeneInfo:
        return "Intergenic","Intergenic", "Intergenic"
    else:
        overlap_ref = GeneInfo[key].find(sidx, eidx)
        overlap_name = [(x.transcriptID, x.gene_name, x.position) for x in overlap_ref if x.position =="gene"]
        if overlap_name:
            return overlap_name[0]
        else:
            overlap_name = [(x.transcriptID, x.gene_name, x.position) for x in overlap_ref]
            if overlap_name:
                return overlap_name[0]
            else:
                return "Intergenic", "Intergenic", "Intergenic"

def main():
    parser = argparse.ArgumentParser(description="anno gene with GRO-cap")
    parser.add_argument("-gene", dest="gene", help="Ghir.gene.bed")
    parser.add_argument("-GWAS", dest="GWAS", help="GWAS.txt")
    parser.add_argument("-o", dest="output", help="output.tsv")
    args = parser.parse_args()
    print(args)
    f_GWAS = args.GWAS
    f_gene = args.gene
    f_out = args.output

    TranscriptInfo = load_transcript_info(f_gene)
    GeneInfo = load_gene(f_gene, flank=1000)
    with open(f_GWAS) as fi, open(f_out, "w") as fo:
        for line in fi:
            if line.startswith("Chrom"):
                header =  line.strip().split("\t") + ["Transcript", "Gene", "Pos", "transPos"]
                fo.write("\t".join(header)+"\n")
                continue
            data = line.strip().split("\t")
            chrom, start, end = data[:3]
            start, end = int(start), int(end)
            overlap_transcript, overlap_gene, overlap_pos = find_gid(chrom, start, end, GeneInfo)
            if overlap_transcript == "Intergenic":
                outline = data + [overlap_transcript, overlap_gene, overlap_pos, "Intergenic"]
            elif overlap_pos != "gene":
                outline = data + [overlap_transcript, overlap_gene, overlap_pos, overlap_pos]
            else:
                trans_info = TranscriptInfo[overlap_transcript]
                trans_pos = trans_info.Region.find(start, end+1)
                # try:
                #     print(trans_pos[0])
                # except:
                #     print(trans_pos, data, overlap_transcript, overlap_gene, overlap_pos)
                outline = data + [overlap_transcript, overlap_gene, overlap_pos, trans_pos[0]]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
