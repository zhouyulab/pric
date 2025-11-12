from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import numpy as np
import datetime
import os,sys, argparse

class Chrimeric(object):
    def __init__(self, chrom, site, strand, block, ReadNum, name, pos):
        self.chrom = chrom
        self.site = int(site)
        self.strand = strand
        self.block = block
        self.ReadNum = int(ReadNum)
        self.name = name
        self.pos = pos

def load_chimeric(f_chimeric):
    ChimericInfo = dict()

    chimeric_n = 0
    with open(f_chimeric) as f:
        for indx, line in enumerate(f):
            if indx == 0:
                continue
            chimeric_n += 1
            chimeric_name = "chimeric%s" % chimeric_n
            data = line.strip().split("\t")
            DonorChrom,DonorCJS, DonorStrand,DonorBlocks,AcceptorChrom,AcceptorCJS, AcceptorStrand, AcceptorBlocks,ReadNum = data
            
            DonorLine = Chrimeric(DonorChrom,DonorCJS, DonorStrand,DonorBlocks, ReadNum, chimeric_name, "Donor")
            if DonorChrom not in ChimericInfo:
                ChimericInfo[DonorChrom] = IntervalTree()
            ChimericInfo[DonorChrom].insert(int(DonorCJS), int(DonorCJS)+1, DonorLine)

            AcceptorLine = Chrimeric(AcceptorChrom,AcceptorCJS, AcceptorStrand, AcceptorBlocks, ReadNum, chimeric_name, "Acceptor")
            if AcceptorChrom not in ChimericInfo:
                ChimericInfo[AcceptorChrom] = IntervalTree()
            ChimericInfo[AcceptorChrom].insert(int(AcceptorCJS), int(AcceptorCJS)+1, AcceptorLine)
    return ChimericInfo


class eQTL(object):
    def __init__(self, line):
        indx, chrom, site, ID, base, group, gene_name = line.strip().split("\t")
        self.chrom = chrom
        self.site = int(site)
        self.ID = ID
        self.group = group
        self.gene_name = gene_name

    def get_eQTL_pos(self, gene_start, gene_end):
        gene_start = int(gene_start)
        gene_end = int(gene_end)
        d1 = self.site - gene_start
        d2 = self.site - gene_end
        diff = np.abs(d1) - np.abs(d2)
        if d1>0 and d2 <0 and diff <0:
            return "up"
        elif d1<=0 and d2 <0:
            return "up"
        elif d1 >0 and d2 <0 and diff>=0:
            return "down"
        elif d1 >0 and d2 >= 0:
            return "down"
        else:
            print(self.start, gene_start, gene_end)

    def get_bins(self, flank, bin, GeneInfo, gene_name):
        gene_chrom, gene_start, gene_end, gene_strand = GeneInfo[gene_name]
        relative_position = self.get_eQTL_pos(gene_start, gene_end)
        if relative_position == "up":
            flank_start = max(self.site-flank, 0)
            flank_end = self.site
        elif relative_position == "down":
            flank_start = self.site
            flank_end = self.site + flank
        flank_bins = np.arange(flank_start, flank_end+bin, step=bin, dtype=int)
        if relative_position == "up":
            flank_bins = flank_bins[::-1]
        return flank_bins
    
def record_Gene(f_gene_bed6):
    GeneInfo = dict()
    with open(f_gene_bed6) as f:
        for line in f:
            data = line.strip().split("\t")
            GeneInfo[data[3]] = [data[0], data[1], data[2], data[-1]]
    return GeneInfo

def record_eQTL(f_eQTL):
    eQTLinfo = dict()
    with open(f_eQTL) as f:
        for line in f:
            data = line.strip().split("\t")
            data = data[1:]
            key = (data[0], data[1])
            if key not in eQTLinfo:
                eQTLinfo[key] = eQTL(line)
    return eQTLinfo

def find_eQTL_bin_chimeric(chrom, bin, gene_name, GeneInfo, ChimericInfo):
    gene_chrom, gene_start, gene_end, gene_strand = GeneInfo[gene_name]
    if (gene_chrom not in ChimericInfo) or (chrom not in ChimericInfo):
        out_bin = [0 for i in range(len(bin)-1)]
        return out_bin
    GeneChimeric = ChimericInfo[gene_chrom].find(int(gene_start), int(gene_end))
    bin_num_li = list()
    for i in range(len(bin)-1):
        bin_start = bin[i]
        bin_end = bin[i+1]
        binChimeric = ChimericInfo[chrom].find(int(bin_start), int(bin_end))
        bin_num = 0
        if GeneChimeric and binChimeric:
            for g_info in GeneChimeric:
                for b_info in binChimeric:
                    if g_info.name == b_info.name and g_info.pos != b_info.pos:
                        assert g_info.ReadNum == b_info.ReadNum
                        bin_num += g_info.ReadNum
        bin_num_li.append(bin_num)
    out_bin = [1 if x>0 else 0 for x in bin_num_li]
    return out_bin

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-e", "--eQTL", type=str, dest="eQTL", metavar="GWAS.integrated.txt", required=True)
    base_group.add_argument("-g", "--gene", type=str, dest="gene", metavar="gene.txt", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_eQTL = args.eQTL
    f_gene = args.gene
    f_out = args.output
    f_chimeric = args.chimeric
    bin = 50
    flank = 2000

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Start....")
    eQTLinfo = record_eQTL(f_eQTL)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Loaded eQTL.")
    GeneInfo = record_Gene(f_gene)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Loaded Gene.")
    ChimericInfo = load_chimeric(f_chimeric)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Loaded Chimeric.")
    with open(f_out, "w") as fo:
        for key in eQTLinfo:
            chrom, site = key
            e_info = eQTLinfo[key]
            gene_name = e_info.gene_name
            eQTL_bins = e_info.get_bins(flank, bin, GeneInfo, gene_name)
            eQTL_bin_bum = find_eQTL_bin_chimeric(chrom, eQTL_bins, gene_name, GeneInfo, ChimericInfo)    
            outline = [chrom, site, gene_name, e_info.group] + eQTL_bin_bum
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()

