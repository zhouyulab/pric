
import sys
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

class SNPline(object):
    def __init__(self, line):
        # print(line)
        chrom, start, end, ref_base, alt_base, trait, Type, SNP, gene_name, region, Reference, Table, refID = line.strip().split("\t")
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.gene_name = gene_name
        self.trait = trait
        self.SNP = SNP
        self.Type = Type
        self.ID = "%s:%s:%s:%s:%s:%s" % (Reference,chrom,start,Type, SNP, trait)

def record_SNP(f_in):
    SNPinfo = dict()
    with open(f_in) as f:
        for line in f:
            if line.startswith("Chrom"):
                continue
            line_info = SNPline(line)
            if line_info.trait == "None":
                continue
            chrom, start, end, ref_base, alt_base, trait, Type, SNP, gene_name, region, Reference, Table, refID = line.strip().split("\t")
            if chrom not in SNPinfo:
                SNPinfo[chrom] = IntervalTree()
            SNPinfo[chrom].insert(int(start), int(end), line_info)

    return SNPinfo

def find_gid(chrom, start, end, GeneInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in GeneInfo:
        return [] 
    overlap_ref = GeneInfo[chrom].find(sidx, eidx)
    overlap_grid = [x.ID for x in overlap_ref]
    return overlap_grid

f_SNP = sys.argv[1]
f_EP = sys.argv[2]
f_out = sys.argv[3]

header = ["Enhancer", "EnhancerPos", "Promoter", "PromoterPos", "SNPloci", "ReadNum", "GeneName", "AtGene", "GeneSymbol", "EnhancerSNP", "PromoterSNP"]
SNPinfo = record_SNP(f_SNP)
with open(f_EP) as f, open(f_out, "w") as f_o:
    f_o.write("\t".join(header)+"\n")
    for line in f:
        data = line.strip().split("\t")
        if data[0] == "Enhancer":
            continue
        enhancer, enhancerPos = data[0], data[1]
        enhancerChrom, enhancerStart, enhancerEnd = enhancerPos.split(":")[0], enhancerPos.split(":")[1].split("-")[0], enhancerPos.split(":")[1].split("-")[1]
        promoter, promoterPos = data[2], data[3]
        promoterChrom, promoterStart, promoterEnd = promoterPos.split(":")[0], promoterPos.split(":")[1].split("-")[0], promoterPos.split(":")[1].split("-")[1]
        enhancerSNP = find_gid(enhancerChrom, enhancerStart, enhancerEnd, SNPinfo)
        promoterSNP = find_gid(promoterChrom, promoterStart, promoterEnd, SNPinfo)

        read_num, gene_name = data[4], data[8]
        AtGene, GeneSymbol = data[9], data[10]
        GeneLabel = "%s(%s)" % (gene_name, GeneSymbol) if  GeneSymbol != "None" else gene_name
        if enhancerSNP and promoterSNP:
            Group = "Both"
        elif enhancerSNP and (not promoterSNP):
            Group = "inEnhancer"
        elif promoterSNP and (not enhancerSNP):
            Group = "inPromoter"
        else:
            Group = "None"

        if Group == "None":
            continue

        enhancer_info = "|".join(enhancerSNP) if enhancerSNP else "None"
        promoter_info = "|".join(promoterSNP) if promoterSNP else "None"
        outline = [enhancer, enhancerPos, promoter, promoterPos, Group, read_num, gene_name, AtGene, GeneLabel, enhancer_info, promoter_info]
        f_o.write("\t".join(outline)+"\n")
