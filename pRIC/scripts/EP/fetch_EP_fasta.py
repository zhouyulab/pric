import pybedtools
import sys, os
import argparse
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

def add_feature(featureInfo, chrom, site, strand, feature):
    key = (chrom, strand)
    site = int(site)
    if key not in featureInfo:
        featureInfo[key].append([feature, int(site)])
        return featureInfo
    is_nearest = False
    for storage_feature, storage_site in featureInfo[key]:
        if abs(storage_site-site)<10:
            is_nearest = True
    if is_nearest:
        return featureInfo
    else:
        featureInfo[key].append([feature, int(site)])
        return featureInfo

def record_feature(f_chimeric):
    PromoterInfo = defaultdict(list)
    EnahncerInfo = defaultdict(list)
    with open(f_chimeric) as fi:
        for line in fi:
            if line.startswith("DonorChrom"):
                continue
            data = line.strip().split("\t")
            DonorChrom, DonorCJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorCJS, AcceptorStrand, AcceptorBlocks, ReadNum,DonorFeature,AcceptorFeature, Group =data
            if Group != "E-P":
                continue
            if DonorFeature.startswith("pro"):
                assert not AcceptorFeature.startswith("pro")
                promoter_chrom, promoter_site, promoter_strand, promoter = DonorChrom, int(DonorCJS), DonorStrand, DonorFeature
                enhancer_chrom, enhancer_site, enhancer_strand, enhancer = AcceptorChrom, int(AcceptorCJS), AcceptorStrand, AcceptorFeature
            else:
                assert AcceptorFeature.startswith("pro")
                assert not DonorFeature.startswith("pro")
                promoter_chrom, promoter_site, promoter_strand, promoter = AcceptorChrom, int(AcceptorCJS), AcceptorStrand, AcceptorFeature
                enhancer_chrom, enhancer_site, enhancer_strand, enhancer = DonorChrom, int(DonorCJS), DonorStrand, DonorFeature
            PromoterInfo = add_feature(PromoterInfo, promoter_chrom, promoter_site, promoter_strand, promoter)
            EnahncerInfo = add_feature(EnahncerInfo, enhancer_chrom, enhancer_site, enhancer_strand, enhancer)
    return PromoterInfo, EnahncerInfo

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fa", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("-pre", "--prefix", type=str, dest="prefix", metavar="out.rprefix", required=True)
    base_group.add_argument("--flank", type=int, dest="flank", metavar=200, required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_genome = args.genome
    f_chimeric = args.chimeric
    f_prefix = args.prefix
    flank = args.flank

    PromoterInfo, EnahncerInfo = record_feature(f_chimeric)
    f_promoter_bed = os.path.join(f_prefix, "promoter.bed")
    f_promoter_fa = os.path.join(f_prefix, "promoter.fa")
    f_enhancer_bed = os.path.join(f_prefix, "enhancer.bed")
    f_enhancer_fa = os.path.join(f_prefix, "enhancer.fa")
    
    with open(f_promoter_bed, "w") as fp:
        for key in PromoterInfo:
            chrom, strand = key
            site_li = PromoterInfo[key]
            for feature, site in site_li:
                flank_start = max(1, site-flank)
                flank_end = site+flank
                outline = [chrom, flank_start, flank_end, feature, site, strand]
                outline = list(map(str, outline))
                fp.write("\t".join(outline)+"\n")
    command1 = "bedtools getfasta -fi %s -bed %s -fo %s -s -name"%(f_genome, f_promoter_bed, f_promoter_fa)
    os.system(command1)

    with open(f_enhancer_bed, "w") as fe:
        for key in EnahncerInfo:
            chrom, strand = key
            site_li = EnahncerInfo[key]
            for feature, site in site_li:
                flank_start = max(1, site-flank)
                flank_end = site+flank
                outline = [chrom, flank_start, flank_end, feature, site, strand]
                outline = list(map(str, outline))
                fe.write("\t".join(outline)+"\n")
    command2 = "bedtools getfasta -fi %s -bed %s -fo %s -s -name"%(f_genome, f_enhancer_bed, f_enhancer_fa)
    os.system(command2)

if __name__ == "__main__":
    main()
