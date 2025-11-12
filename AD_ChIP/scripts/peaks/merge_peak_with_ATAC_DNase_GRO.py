import pybedtools
import sys
import argparse
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

def write_promoter(f_in, f_out):
    GenePromoter = defaultdict(list)
    PromoterInfo = dict()
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split("\t")
            GenePromoter[data[-2]].append(data)
    with open(f_out, "w") as fo:
        for gene in GenePromoter:
            gene_info = GenePromoter[gene]
            start_li = [int(x[1]) for x in gene_info]
            end_li = [int(x[2]) for x in gene_info]
            sample_li = list()
            for x in gene_info:
                sample_li.extend(x[3].split(","))
            sample_li =  sorted(list(set(sample_li)))
            outline = [gene_info[0][0], min(start_li), max(end_li), "gro%s" % gene, ",".join(sample_li), "."]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

            key = gene_info[0][0]
            if key not in PromoterInfo:
                PromoterInfo[key] = IntervalTree()
            PromoterInfo[key].insert(min(start_li), max(end_li), Bedline("\t".join(outline)+"\n"))
    return PromoterInfo

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

def load_peak_info(f_gro):
    GROInfo = dict()
    with open(f_gro) as f:
        for line in f:
            groline = Bedline(line)
            chrom, start, end, name, count ,strand = line.strip().split("\t")
            start, end = int(start), int(end)
            key = chrom
            if key not in GROInfo:
                GROInfo[key] = IntervalTree()
            GROInfo[key].insert(start, end, groline)
    return GROInfo


def record_ChIP(f_ChIP):
    ChIP_info = dict()
    with open(f_ChIP) as fi:
        for line in fi:
            data = line.strip().split("\t")
            ChIP_info[data[4]] = data
    return ChIP_info

def recorf_merge_info(Peak_merge):
    merge_Info = defaultdict(list)
    for x in Peak_merge:
        key = x[4]
        merge_Info[key].append(x)
    return merge_Info

def reshape_pos(start_li, end_li, pos_info):
    for pos in pos_info:
        start = int(pos[7])
        end = int(pos[8])
        start_li.append(int(start))
        end_li.append(int(end))
    return start_li, end_li

def find_gid(chrom, start, end, GROInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in GROInfo:
        return [] 
    overlap_ref = GROInfo[chrom].find(sidx, eidx)
    if not overlap_ref:
        return [] 
    gro_strand = [x.strand for x in overlap_ref]
    if "+" in gro_strand and "-" in gro_strand:
        gro_group = "both"
    elif "+" in gro_strand:
        gro_group = "forward"
    elif "-" in gro_strand:
        gro_group = "reverse"
    else:
        return "."
    return gro_group

def record_all_info(f_ChIP, f_DNase, f_ATAC):
    PeakInfo = list()
    with open(f_ChIP) as fi:
        for line in fi:
            data = line.strip().split("\t")
            chrom, start, end, sample, name, strand = data
            chip_line = [chrom, start, end, "ChIP", name, strand]
            PeakInfo.append(chip_line)
    with open(f_ChIP) as fi:
        for line in fi:
            data = line.strip().split("\t")
            chrom, start, end, sample, name, strand = data
            chip_line = [chrom, start, end, "ChIP", name, strand]
            PeakInfo.append(chip_line)


    return ChIP_info


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("--DNase", type=str, dest="DNase", metavar="DNase.peak", required=True)
    base_group.add_argument("--ATAC", type=str, dest="ATAC", metavar="ATAC.peak", required=True)
    base_group.add_argument("--H3K4me4", type=str, dest="H3K4me4", metavar="H3K4me4.peak", required=True)
    base_group.add_argument("--H3K27ac", type=str, dest="H3K27ac", metavar="H3K27ac.peak", required=True)
    base_group.add_argument("--GRO", type=str, dest="GRO", metavar="GRO.peak", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.bed", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_DNase = args.DNase
    f_ATAC = args.ATAC
    f_H3K27ac = args.H3K27ac
    f_H3K4me4 = args.H3K4me4
    f_GRO = args.GRO
    f_enahncer = args.enhancer
    f_promoter = args.promoter

    # PromoterInfo = write_promoter(f_H3K4me4, f_promoter) 

    Peak_Info = record_ChIP(f_H3K27ac)
    ChIP_Peak = pybedtools.BedTool(f_H3K27ac)
    ATAC_Peak = pybedtools.BedTool(f_ATAC)
    DNase_Peak = pybedtools.BedTool(f_DNase)
    H3K4me4_Peak = pybedtools.BedTool(f_H3K4me4)
    GROInfo = load_peak_info(f_GRO)

    ChIP_ATAC_merge = recorf_merge_info(ChIP_Peak.intersect(ATAC_Peak, wo=True))
    ChIP_DNase_merge = recorf_merge_info(ChIP_Peak.intersect(DNase_Peak, wo=True))
    ChIP_H3K4me4_merge = recorf_merge_info(ChIP_Peak.intersect(H3K4me4_Peak, wo=True))

    print(len(ChIP_ATAC_merge), len(ChIP_DNase_merge), len(ChIP_H3K4me4_merge))
    Peaks = list()
    for key in Peak_Info:
        atac_info = ChIP_ATAC_merge[key]
        dnase_info = ChIP_DNase_merge[key]
        chip_info = Peak_Info[key]
        promoter_info = ChIP_H3K4me4_merge[key]
        chrom, start, end = chip_info[0], int(chip_info[1]), int(chip_info[2])
        sample_li = chip_info[3]

        if promoter_info:
            continue
        if not dnase_info and not atac_info:
            continue
        
        start_li, end_li = [start], [end]
        if dnase_info:
            start_li, end_li = reshape_pos(start_li, end_li, dnase_info)
        if atac_info:
            start_li, end_li = reshape_pos(start_li, end_li, atac_info)
        
        enhancer_start = min(start_li)
        enahncer_end = max(end_li)

        gro_group = find_gid(chrom, enhancer_start, enahncer_end, GROInfo)
        promoter_info = find_gid(chrom, enhancer_start, enahncer_end, PromoterInfo)
        if gro_group and not promoter_info:
            peak_line = [chrom, enhancer_start, enahncer_end, sample_li, gro_group, "."]
            Peaks.append(peak_line)
    Peaks = pybedtools.BedTool(Peaks)
    Peak_sort = Peaks.sort().merge( d=200, c='4,5', o="distinct")
    
    n = 0
    with open(f_enahncer, "w") as fo:
        for x in Peak_sort:
            n += 1
            enhancer_name = "enhancer%s" % n
            enhancer_sample = ",".join(list(set(x[3].split(","))))
            if "both" in x[4]:
                enhancer_strand = "both"
            elif "forward" in x[4] and "reverse" in x[4]:
                enhancer_strand = "both"
            else:
                enhancer_strand = list(set(x[4].split(",")))
                assert len(enhancer_strand) == 1
                enhancer_strand = enhancer_strand[0]
            outline = [x[0], x[1], x[2], enhancer_name, enhancer_sample, ".", enhancer_strand]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
