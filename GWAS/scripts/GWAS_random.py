from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import numpy as np
import datetime
import random
import os,sys, argparse
random.seed(123)

class GWAS(object):
    def __init__(self, line):
        Chrom,start, end, ref_base, alt_base, trait, Type, SNP, gene_name, region, reference, table, refID = line.strip().split("\t")
        self.chrom = Chrom
        self.site = int(start)
        self.SNP = SNP
        self.refBase = ref_base
        self.altBase = alt_base
        self.Trait = trait
        self.Gene = gene_name
        self.Source = reference

def record_GWAS(f_GWAS):
    GWASinfo = dict()
    with open(f_GWAS) as f:
        for line in f:
            if line.startswith("Chrom"):
                continue
            data = line.strip().split("\t")
            Chrom,start, end, ref_base, alt_base, trait, Type, SNP, gene_name, region, reference, table, refID = data

            GWASLine = GWAS(line)
            if Chrom not in GWASinfo:
                GWASinfo[Chrom] = IntervalTree()
            GWASinfo[Chrom].insert(int(start), int(end)+1, GWASLine)
    return GWASinfo

class Junction(object):
    def __init__(self, Chrom, Site):
        self.chrom = Chrom
        self.site = int(Site)
        self.bins = None

    def get_flank_bin(self, flank, bin):
        flank_start = max(self.site-flank, 0)
        flank_end = self.site+flank 
        flank_bins = np.arange(flank_start, flank_end+bin, step=bin, dtype=int)
        self.bins = flank_bins
        return flank_bins

    def get_GWAS_info(self, flank, bin, GWASinfo):
        if self.bins is None:
            self.get_flank_bin(flank, bin)

        bin_num_li = list()
        for i in range(len(self.bins)-1):
            bin_start = self.bins[i]
            bin_end = self.bins[i+1]
            if self.chrom not in GWASinfo:
                bin_num = 0
            else:
                binGWAS = GWASinfo[self.chrom].find(int(bin_start), int(bin_end))
                if binGWAS:
                    bin_num = 1
                else:
                    bin_num = 0
            bin_num_li.append(bin_num)
        return bin_num_li
        

def record_junction(f_chimeric):
    JunctionInfo = dict()
    with open(f_chimeric) as f:
        for line in f:
            if line.startswith("DonorChrom"):
                continue
            data = line.strip().split("\t")
            DonorChrom, DonorJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorJS, AcceptorStrand, AcceptorBlocks, ReadNum = data  
            if  DonorChrom not in JunctionInfo:
                JunctionInfo[DonorChrom] = IntervalTree()
            JunctionInfo[DonorChrom].insert(max(int(DonorJS)-2000, 0), int(DonorJS)+2000, "Donor")

            if AcceptorChrom not in JunctionInfo:
                JunctionInfo[AcceptorChrom] = IntervalTree()
            JunctionInfo[AcceptorChrom].insert(max(int(AcceptorJS)-2000, 0), int(AcceptorJS)+2000, "Acceptor")
    return JunctionInfo

def record_chrom_size(f_chrom_size):
    ChromSizeDict = dict()
    with open(f_chrom_size) as f:
        for line in f:
            data = line.strip().split("\t")
            if data[0].startswith("Scaffold"):
                continue
            ChromSizeDict[data[0]] = int(data[1])
    return  ChromSizeDict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("-g", "--GWAS", type=str, dest="GWAS", metavar="GWAS.txt", required=True)
    base_group.add_argument("-chrom", "--chrom-size", type=str, dest="chrom_size", metavar="chrom_size.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_chimeric = args.chimeric
    f_GWAS = args.GWAS
    f_out = args.output
    f_chromSize = args.chrom_size

    JunctionInfo = record_junction(f_chimeric)
    GWASinfo = record_GWAS(f_GWAS)
    ChromSizeDict = record_chrom_size(f_chromSize)
    Number = 7903322
    flank=1000
    bin = 50

    out_li  = [0 for i in range(40)]
    total_read = 0
    for i in range(Number):
        isOverlap = True
        while isOverlap:
            random_chrom = random.sample(ChromSizeDict.keys(), 1)[0]
            random_site = random.randrange(1,ChromSizeDict[random_chrom])
            res_overlap = JunctionInfo[random_chrom].find(random_site, random_site+1)
            if not res_overlap:
                isOverlap = False
        random_info = Junction(random_chrom, random_site)
        random_bin = random_info.get_GWAS_info(flank, bin, GWASinfo)
        if len(random_bin) != 40:
            continue
        total_read += 1
        out_li = [out_li[i]+random_bin[i] for i in range(40)]

    header = ["Index", "Num", "Group"]
    with open(f_out, "w") as fo:
        fo.write("\t".join(header)+"\n")
        total_line = [str(total_read), str(total_read), "total"]
        fo.write("\t".join(total_line)+"\n")
        for indx, num in enumerate(out_li):
            outline = [indx, num, "random"]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")
if __name__ == "__main__":
    main()


    # f_chimeric = "../../RIC-seq/results/chimeric/merge/merge_ChimeReadCnt.tsv"
    # f_GWAS = "../results/GWAS.mapping/GWAS.integrated.txt"
    # f_chrom_size = "../../genome/Ghir_genome/Ghirsutum.chrom.sizes"