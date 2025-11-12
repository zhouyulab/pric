#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

class Bedline(object):
    def __init__(self, line):
        # print(line)
        chrom, start, end, name, source, strand = line.strip().split("\t")
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand


def load_gene_bed6(genebed):
    GeneInfo = dict()
    with open(genebed) as f:
        for line in f:
            bedline = Bedline(line)
            line = line.strip().split("\t")
            chrom, start, end, name, source, strand = line
            gene_type, gene_name = name.split("|")
            if gene_type in ["TE"]:
                continue
            key = (chrom, strand)
            if key not in GeneInfo:
                GeneInfo[key] = IntervalTree()
            GeneInfo[key].insert(int(start), int(end), bedline)
    return GeneInfo

def find_gid(chrom, site, strand, GeneInfo):
    sidx = int(site)
    eidx = int(site)+1
    if (chrom, strand) not in GeneInfo:
        return [] 
    overlap_ref = GeneInfo[(chrom, strand)].find(sidx, eidx)
    overlap_gids = [x.name for x in overlap_ref]
    return overlap_gids


def gene_interaction(inputtsv, GeneInfo, outputtsv):
    with open(outputtsv, 'w') as fo, open(inputtsv) as f_in:
        for indx, line in enumerate(f_in):
            if indx == 0:
                header = line.strip().split("\t") + ["DonorGene", "AccetperGene", "Group"]
                fo.write("\t".join(header)+"\n")
                continue
            
            data = line.rstrip("\n").split("\t")        
            DonorChrom, DonorJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorJS, AcceptorStrand, AcceptorBlocks, ReadNum = data   
            DonorGids = find_gid(DonorChrom, DonorJS, DonorStrand, GeneInfo)
            AcceptorGids = find_gid(AcceptorChrom, AcceptorJS, AcceptorStrand, GeneInfo)
            if DonorGids and AcceptorGids:
                DonorGene = ";".join(DonorGids)
                AcceptorGene = ";".join(AcceptorGids)
                if list(set(DonorGids) & set(AcceptorGids)):
                    group = "intra-gene"
                else:
                    group = "inter-gene"
            elif DonorGids and not AcceptorGids:
                DonorGene = ";".join(DonorGids)
                AcceptorGene = "None"
                group = "gene-None"
            elif not DonorGids and AcceptorGids:
                DonorGene = "None"
                AcceptorGene = ";".join(AcceptorGids)
                group = "gene-None"
            else:
                DonorGene = "None"
                AcceptorGene = "None"
                group = "None-None"
            outline = data + [DonorGene, AcceptorGene, group]
            fo.write("\t".join(outline)+"\n")

def statis_group_num(f_in, f_out):
    GroupStatis = defaultdict(int)
    with open(f_in) as f:
        for line in f:
            if line.startswith("DonorChrom"):
                continue
            data =line.strip().split("\t")
            group, readnum = data[-1], int(data[-4])
            GroupStatis[group] += readnum
    with open(f_out, "w") as f_o:
        for group in GroupStatis:
            outline = [group, str(GroupStatis[group])]
            f_o.write("\t".join(outline)+"\n")

def main():
    parser = argparse.ArgumentParser(description="Fetch gene interactions")
    parser.add_argument("-i", dest="input", help="ChimericRead.tsv")
    parser.add_argument("-b", dest="ref", help="gene.bed6")
    parser.add_argument("-o", dest="output", help="ChimJS.gene_interaction.cnt.tsv")
    parser.add_argument("-og", dest="outstaris", help="ChimJS.Group.tsv")
    args = parser.parse_args()
    f_in = args.input
    f_genebed = args.ref
    f_out = args.output
    f_outstatis = args.outstaris
    GeneInfo = load_gene_bed6(f_genebed)
    gene_interaction(f_in, GeneInfo, f_out)
    statis_group_num(f_out,f_outstatis)

if __name__ == "__main__":
    main()
