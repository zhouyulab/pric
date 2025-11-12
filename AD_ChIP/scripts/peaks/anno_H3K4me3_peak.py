import pybedtools
import sys
import argparse
from collections import defaultdict
from bx.intervals.intersection import IntervalTree
from interval import Interval

class Bedline(object):
    def __init__(self, line):
        # print(line)
        chrom, start, end, name, source, strand = line.strip().split("\t")
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand

def load_gene_bed6(genebed, flank):
    GeneInfo = dict()
    with open(genebed) as f:
        for line in f:
            bedline = Bedline(line)
            line = line.strip().split("\t")
            chrom, start, end, name, source, strand = line
            start, end = int(start), int(end)
            if strand == "+":
                flank_start = max(1, start-flank)
                flank_end = start + flank
            else:
                flank_start = max(1, end-flank)
                flank_end = end + flank
            key = chrom
            if key not in GeneInfo:
                GeneInfo[key] = IntervalTree()
            GeneInfo[key].insert(int(flank_start), int(flank_end), bedline)
    return GeneInfo

def get_overlap_length(I1, I2):
    if I2 in I1:
        length = I2.upper_bound - I2.lower_bound
    elif I1 in I2:
        length = I2.upper_bound - I1.lower_bound
    else:
        I3 = I1 & I2
        length = I3.upper_bound - I3.lower_bound
    # assert length > 0
    return length


def find_gid(chrom, start, end, GeneInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in GeneInfo:
        return [] 
    overlap_ref = GeneInfo[chrom].find(sidx, eidx)
    overlap_length = list()
    for x in overlap_ref:
        I1 = Interval(x.start, x.end)
        I2 = Interval(sidx, eidx)
        if I1.overlaps(I2):
            length =  get_overlap_length(I1, I2)
            overlap_length.append((x.name, length))
    if overlap_length:
        overlap_length.sort(key=lambda x: int(x[1]), reverse=True)
        return overlap_length[0][0]
    else:
        return []
def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--peak", type=str, dest="peak", metavar="Hak4me3_merge.peak", required=True)
    base_group.add_argument("-b", "--bed", type=str, dest="bed", metavar="Gh.gene.bed6", required=True)
    base_group.add_argument("-f", "--flank", type=int, dest="flank", metavar=2000, required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="outdir", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_peak = args.peak
    f_bed = args.bed
    flank = args.flank
    f_out = args.output

    BedInfo = load_gene_bed6(f_bed, flank)
    n = 0
    with open(f_peak) as fi, open(f_out, "w") as fo:
        for line in fi:
            chrom, start, end, sample, name, strand = line.strip().split("\t")
            gene_name = find_gid(chrom, start, end, BedInfo)
            if gene_name:
                n = n+1
                promoter_name = "promoter%s" % n 
                outline = [chrom, start, end, sample, name, strand, gene_name, promoter_name]
                fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()