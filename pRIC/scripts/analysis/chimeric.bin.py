import numpy as np
import argparse
import os
import sys
from collections import defaultdict
from bx.intervals.intersection import IntervalTree


def record_chrom_size(f_chrom):
    ChromSize = dict()
    with open(f_chrom) as f:
        for line in f:
            data = line.strip().split("\t")
            chrom, size = data[0], int(data[1])
            if "Scaffold" in chrom:
                continue
            ChromSize[chrom] = size
    return ChromSize

class Bedline(object):
    def __init__(self, chrom, start, end, name):
        # print(line)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name


def record_bin_info(ChromSize, resolution, outid):
    sort_chrims = sorted(ChromSize)
    bin_id = 0
    GeneInfo = dict()
    IDinfo = dict()
    header = ["Chrom", "Start", "End", "ID"]
    with open(outid, "w") as f_o:
        f_o.write("\t".join(header)+"\n")
        for chrom in sort_chrims:
            size = ChromSize[chrom]
            size_li = np.arange(0, size, dtype=int, step=resolution)
            GeneInfo[chrom] = IntervalTree()
            for i in range(len(size_li)-1):
                bin_id = bin_id + 1
                outline = [chrom, size_li[i], size_li[i+1], bin_id]
                f_o.write("\t".join(list(map(str,outline)))+"\n")
                IDinfo[bin_id] = outline
                bed_line = Bedline(chrom, size_li[i], size_li[i+1], bin_id)
                GeneInfo[chrom].insert(int(size_li[i]), int(size_li[i+1]), bed_line)
                
            bin_id = bin_id + 1     
            outline = [chrom, size_li[-1], size, bin_id]
            f_o.write("\t".join(list(map(str,outline)))+"\n")
            IDinfo[bin_id] = outline
            bed_line = Bedline(chrom, size_li[-1], size, bin_id)
            GeneInfo[chrom].insert(int(size_li[-1]), size, bed_line)
    return GeneInfo, IDinfo

def find_gid(chrom, site, strand, GeneInfo):
    sidx = int(site)
    eidx = int(site)+1
    if chrom not in GeneInfo:
        return [] 
    overlap_ref = GeneInfo[chrom].find(sidx, eidx)
    overlap_gids = [x.name for x in overlap_ref]
    return overlap_gids

def record_read_info(inputtsv, GeneInfo):
    ReadsInfo = defaultdict(int)
    with  open(inputtsv) as f_in:
        for indx, line in enumerate(f_in):
            if indx == 0:
                continue
            data = line.rstrip("\n").split("\t")        
            DonorChrom, DonorJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorJS, AcceptorStrand, AcceptorBlocks, ReadNum = data   
            DonorGids = find_gid(DonorChrom, DonorJS, DonorStrand, GeneInfo)
            AcceptorGids = find_gid(AcceptorChrom, AcceptorJS, AcceptorStrand, GeneInfo)
            read_num = int(ReadNum)
            for d_grid in DonorGids:
                for a_grid in AcceptorGids:
                    if d_grid < a_grid:
                        key = (d_grid, a_grid)
                    else:
                        key = (a_grid, d_grid)
                    ReadsInfo[key] += read_num
    return ReadsInfo

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-c", "--chrom-size", type=str, dest="chrom_size", metavar="chrom.size.txt", required=True)
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="ChimeReadCnt.tsv", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("--out-id", "--output-id", type=str, dest="outputID", metavar="output.txt", required=True)
    base_group.add_argument("-r", "--resolution", type=int, dest="resolution", metavar=1000000, required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_chrom = args.chrom_size
    inputtsv = args.input
    f_out = args.output
    f_outid = args.outputID
    resolution = args.resolution

    ChromSize = record_chrom_size(f_chrom)
    GeneInfo, IDinfo = record_bin_info(ChromSize, resolution, f_outid)
    ReadsInfo = record_read_info(inputtsv, GeneInfo)
    print(len(ReadsInfo))
    header = ["ID1", "Chrom1", "Start1", "End1", "ID2", "Chrom2", "Start2", "End2", "ReadNum"]
    with open(f_out, "w") as f_o:
        f_o.write("\t".join(header)+"\n")
        for d_grid, a_grid in ReadsInfo:
            chrom1, start1, end1, id1 = IDinfo[d_grid]
            chrom2, start2, end2, id2 = IDinfo[a_grid]
            read_num = ReadsInfo[(d_grid, a_grid)]
            outline = [d_grid, chrom1, start1, end1, a_grid, chrom2, start2, end2, read_num]
            outline = list(map(str, outline))
            f_o.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()