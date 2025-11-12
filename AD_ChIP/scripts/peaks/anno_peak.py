import pybedtools
import sys
import argparse
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

class Bedline(object):
    def __init__(self, line):
        # print(line)
        chrom, start, end, name, count ,strand = line.strip().split("\t")
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand
        self.pos = "%s:%s-%s:%s" % (chrom, start, end, name)

def load_bed_info(f_in):
    GeneInfo = dict()
    with open(f_in) as f:
        for line in f:
            proline = Bedline(line)
            chrom, start, end, name, one ,strand = line.strip().split("\t")
            start, end = int(start), int(end)
            key = chrom
            if key not in GeneInfo:
                GeneInfo[key] = IntervalTree()
            GeneInfo[key].insert(start, end, proline)
    return GeneInfo

def find_gid(chrom, start, end, GeneInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in GeneInfo:
        return "None"
    overlap_ref = GeneInfo[chrom].find(sidx, eidx)
    if overlap_ref:
        overlap_gene = [x.name for x in overlap_ref]
        return ";".join(overlap_gene)
    else:
        return "Intergenic"


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="merge_all.gene.bed6", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_bed = args.ref
    f_enahncer = args.enhancer
    f_output = args.output

    GeneInfo = load_bed_info(f_bed)
    with open(f_enahncer) as fi, open(f_output, "w") as fo:
        for line in fi:
            chrom, start,end, enhancer, enhancer_sample, strand, enhancer_group = line.strip().split("\t")
            start,end = int(start), int(end)
            gene_anno = find_gid(chrom, start, end, GeneInfo)
            outline = line.strip().split("\t") + [gene_anno]
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
