import pybedtools
import sys
import argparse
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

class Bedline(object):
    def __init__(self, chrom, start, end, name, sample, strand, feature):
        # print(line)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand
        self.sample = sample
        self.feature = feature

    def get_info(self):
        return [self.chrom, self.start, self.end, self.name, self.feature, self.strand]

def load_feature_bed6(f_promoter, f_enhancer):
    FeatureInfo = dict()
    with open(f_enhancer) as f:
        for line in f:
            line = line.strip().split("\t")
            chrom, start, end, name, sample, strand, pro = line
            bedline = Bedline(chrom, start, end, name, sample, strand, "enhancer")
            start, end = int(start), int(end)
            key = chrom
            if key not in FeatureInfo:
                FeatureInfo[key] = IntervalTree()
            FeatureInfo[key].insert(start, end, bedline)

    with open(f_promoter) as f:
        for line in f:
            line = line.strip().split("\t")
            chrom, start, end, name, sample, strand = line
            bedline = Bedline(chrom, start, end, name, sample, strand, "promoter")
            start, end = int(start), int(end)
            key = chrom
            if key not in FeatureInfo:
                FeatureInfo[key] = IntervalTree()
            FeatureInfo[key].insert(start, end, bedline)
    return FeatureInfo

def find_gid(chrom, start, end, GeneInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in GeneInfo:
        return [] 
    overlap_ref = GeneInfo[chrom].find(sidx, eidx)
    return overlap_ref


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.peak", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.peak", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("-n", "--network", type=str, dest="network", metavar="EP.network", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_promoter = args.promoter
    f_enhancer = args.enhancer
    f_output = args.output
    f_chimeric = args.chimeric
    f_network = args.network

    FeatureInfo = load_feature_bed6(f_promoter, f_enhancer)
    EPcount = defaultdict(int)
    with open(f_chimeric) as fi, open(f_output, "w") as fo:
        for line in fi:
            if line.startswith("DonorChrom"):
                header = line.strip().split("\t") + ["DonorFeature", "AcceptorFeature", "Group"]
                fo.write("\t".join(header)+"\n")
                continue
            data = line.strip().split("\t")
            DonorChrom, DonorCJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorCJS, AcceptorStrand, AcceptorBlocks, ReadNum =data
            DonorCJS, AcceptorCJS = int(DonorCJS), int(AcceptorCJS)
            DonorGid = find_gid(DonorChrom, DonorCJS, DonorCJS+1, FeatureInfo)
            AcceptorGid = find_gid(AcceptorChrom, AcceptorCJS, AcceptorCJS+1, FeatureInfo)
            if DonorGid and AcceptorGid:
                for d_gid in DonorGid:
                    for a_gid in AcceptorGid:
                        if d_gid.feature == "promoter" and a_gid.feature == "promoter":
                            group = "P-P"
                        elif d_gid.feature == "enhancer" and a_gid.feature == "enhancer":
                            group = "E-E"
                        else:
                            group = "E-P"
                        d_gid_info = d_gid.get_info()
                        a_gid_info = a_gid.get_info()
                        key = tuple(d_gid_info+a_gid_info) if a_gid.name < d_gid.name else tuple(a_gid_info+d_gid_info)
                        EPcount[key] += int(ReadNum)
                        outline = data + [d_gid.name, a_gid.name, group]
                        fo.write("\t".join(outline)+"\n")

    with open(f_network, "w") as fn:
        for key in EPcount:
            reads = EPcount[key]
            outline = list(key) + [reads]
            outline = list(map(str, outline))
            fn.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
