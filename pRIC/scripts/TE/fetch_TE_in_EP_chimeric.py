from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import numpy as np
import datetime
import os,sys, argparse


class TEInfo(object):
    def __init__(self, line):
        data = line.strip().split("\t")
        chrom, start, end, repeat, family = data[4], data[5], data[6], data[9], data[10]
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.repeat = repeat
        self.TE = family
        self.info = "%s::%s:%s-%s" % (self.TE, self.chrom, start, end)
    
    def get_TE_overlap_ratio(self, chrom, start, end, cutoff = 0.5):
        assert self.chrom == chrom
        assert self.start < self.end
        assert start < end
        if start <= self.start and end >= self.end:
            return True
        elif self.start <= start and self.end >= end:
            return True
        elif start <= self.start:
            overlap_length = end - self.start
            ratio1 = overlap_length / (end-start)
            ratio2 = overlap_length / (self.end-self.start)
            if max(ratio1, ratio2) > cutoff:
                return True
            else:
                return False
        else:
            overlap_length = self.end - start
            ratio1 = overlap_length / (end-start)
            ratio2 = overlap_length / (self.end-self.start)
            if max(ratio1, ratio2) > cutoff:
                return True
            else:
                return False

def record_TE(f_bed):
    TEbindInfo = dict()
    with open(f_bed) as f:
        for line in f:
            data = line.strip().split("\t")
            if data[0] == "score":
                continue
            chrom, start, end, repeat, family = data[4], data[5], data[6], data[9], data[10]
            TE_line = TEInfo(line)
            if chrom not in TEbindInfo:
                TEbindInfo[chrom] = IntervalTree()
            TEbindInfo[chrom].insert(int(start), int(end), TE_line)
    return TEbindInfo


def record_chimeric(f_in):
    ChimericInfo = defaultdict(int)
    total_reads = 0
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split("\t")
            if data[0] == "DonorChrom":
                continue
            DonorChrom,DonorCJS,DonorStrand,DonorBlocks,AcceptorChrom,AcceptorCJS, AcceptorStrand,AcceptorBlocks,ReadNum,DonorFeature,AcceptorFeature, Group = data
            ReadNum = int(ReadNum)
            if Group != "E-P":
                continue
            if DonorFeature.startswith("pro"):
                assert not AcceptorFeature.startswith("pro")
                key = (AcceptorFeature, AcceptorChrom, AcceptorCJS, DonorFeature, DonorChrom, DonorCJS)
            else:
                assert AcceptorFeature.startswith("pro")
                key = (DonorFeature, DonorChrom, DonorCJS, AcceptorFeature, AcceptorChrom, AcceptorCJS)
            total_reads += ReadNum
            ChimericInfo[key] += ReadNum
    return ChimericInfo, total_reads


def find_TE_info(Chrom, Site, flank, TEbindInfo):
    chrom = Chrom
    start = max(1, int(Site)-flank)
    end = int(Site)+flank
    if chrom not in TEbindInfo:
        return []
    TE_pos = TEbindInfo[chrom].find(int(start), int(end))
    TE_ref = [x.TE for x in TE_pos]
    TE_overlap = [x.get_TE_overlap_ratio(chrom, start, end, cutoff = 0.5) for x in TE_pos]
    TE_filter = [TE_ref[i] for i in range(len(TE_ref)) if TE_overlap[i] ]
    return TE_filter


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-TE", "--transposon", type=str, dest="transposon", metavar="genome.fa.out.txt", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("--flank", type=int, dest="flank", metavar=200, required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("-s", "--statis", type=str, dest="statis", metavar="statis.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_chimeric = args.chimeric
    f_TE = args.transposon
    f_out = args.output
    flank = args.flank
    f_statis = args.statis
    
    TEbindInfo = record_TE(f_TE)
    ChimericInfo, total_reads = record_chimeric(f_chimeric)

    TE_contain_enhancer = defaultdict(int)
    TE_contain_promoter = defaultdict(int)
    header = ["Enhancer", "EnhancerChrom", "EnhancerSite", "Promoter", "PromoterChrom", "PromoterSite", "Reads", "Total", "EnhancerTE", "PromoterTE", "TEGroup"]
    with open(f_out, "w") as fo:
        fo.write("\t".join(header)+"\n")
        for key in ChimericInfo:
            reads = ChimericInfo[key]
            enahncer, enhancer_chrom, enhancer_site, promoter, promoter_chrom, promoter_site = key
            EnhancerTE = find_TE_info(enhancer_chrom, enhancer_site, flank, TEbindInfo)
            PromoterTE = find_TE_info(promoter_chrom, promoter_site, flank, TEbindInfo)
            if EnhancerTE and PromoterTE:
                TE_group = "EP-pair"
            elif EnhancerTE and not PromoterTE:
                TE_group = "enhancer"
            elif not EnhancerTE and PromoterTE:
                TE_group = "promoter"
            else:
                TE_group = "non-TE"
            outline = [enahncer, enhancer_chrom, enhancer_site, promoter, promoter_chrom, promoter_site, str(reads), str(total_reads), "|".join(EnhancerTE), "|".join(PromoterTE), TE_group]
            fo.write("\t".join(outline)+"\n")
        

            for te in list(set(EnhancerTE)):
                TE_contain_enhancer[te] += reads
            for te in list(set(PromoterTE)):
                TE_contain_promoter[te] += reads

    with open(f_statis, "w") as fo:
        fo.write("Feature\tTransposable\tNum\tTotal\tPercent\n")
        for te in TE_contain_enhancer:
            te_num = TE_contain_enhancer[te]
            percent = te_num * 100 /total_reads
            outline = ["enhancer", te, te_num, total_reads, round(percent, 4)]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")
        for te in TE_contain_promoter:
            te_num = TE_contain_promoter[te]
            percent = te_num * 100 /total_reads
            outline = ["promoter", te, te_num, total_reads, round(percent, 4)]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")


if __name__ == "__main__":
    main()
