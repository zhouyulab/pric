import pybedtools
import sys,os
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
            outline = [gene_info[0][0], min(start_li), max(end_li), "pro%s" % gene, ",".join(sample_li), "."]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

            key = gene_info[0][0]
            if key not in PromoterInfo:
                PromoterInfo[key] = IntervalTree()
            PromoterInfo[key].insert(min(start_li), max(end_li), Bedline("\t".join(outline)+"\n"))
    return PromoterInfo

def merge_all_promoter(f_bed, PromoterInfo, flank=2000):
    with open(f_bed) as f:
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
            if key not in PromoterInfo:
                PromoterInfo[key] = IntervalTree()
            PromoterInfo[key].insert(int(flank_start), int(flank_end), bedline)
    return PromoterInfo


def merge_ChIP_peak(f_ChIP, PromoterInfo):
    PeakInfo = list()
    with open(f_ChIP) as fi:
        for line in fi:
            data = line.strip().split("\t")
            chrom, start, end, name, length, strand = data[:6]
            promoter_info = find_gid(chrom, start, end, PromoterInfo)
            if promoter_info:
                continue
            chip_line = [chrom, start, end, "ChIP", name, strand]
            PeakInfo.append(chip_line)
    Peaks = pybedtools.BedTool(PeakInfo)
    Peak_sort = Peaks.sort().merge(d=1000, c='4,5', o="distinct")
    return Peak_sort

def load_GRO_info(f_gro):
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


def record_ChIP_DNase(ChIPPeak, f_DNase):
    PeakInfo = list()
    for x in ChIPPeak:
        chrom, start, end, chip, sample = x
        strand = "."
        chip_line = [chrom, start, end, "ChIP", sample, strand]
        PeakInfo.append(chip_line)
        
    with open(f_DNase) as fi:
        for line in fi:
            data = line.strip().split("\t")
            chrom, start, end, sample, name, strand = data
            dnase_line = [chrom, start, end, "DNase", name, strand]
            PeakInfo.append(dnase_line)

    Peaks = pybedtools.BedTool(PeakInfo)
    Peak_sort = Peaks.sort().merge(d=0, c='4,5', o="distinct")

    return Peak_sort

def class_peak_group(merge_info):
    if ("ChIP" in merge_info) and (("DNase" in merge_info) or ("ATAC" in merge_info)):
        peak_group = "H3K27ac_Chromatin"
    elif ("ChIP" in merge_info) and ("DNase" not in merge_info) and ("ATAC" not in merge_info):
        peak_group = "H3K27ac"
    elif ("ChIP" not in merge_info) and (("DNase" in merge_info) or ("ATAC" in merge_info)):
        peak_group = "Chromatin"
    else:
        peak_group = "None"
    return peak_group

def group_peak(PeakInfo, PromoterInfo):
    merge_Peak = list()
    H3K27ac_Peak = list()
    chromatin_Peak = list()
    PeakInfo2 = PeakInfo.sort().merge(d=500, c='4,5', o="distinct")
    for x in PeakInfo2:
        chrom, start, end, group, sample = x
        start, end = int(start), int(end)
        promoter_info = find_gid(chrom, start, end, PromoterInfo)
        peak_group = class_peak_group(group)
        if promoter_info:
            continue
        peak_line = [chrom, start, end, group, sample]
        if peak_group == "H3K27ac_Chromatin":
            merge_Peak.append(peak_line)
        elif peak_group == "H3K27ac":
            H3K27ac_Peak.append(peak_line)
        elif peak_group == "Chromatin":
            chromatin_Peak.append(peak_line)
    return pybedtools.BedTool(merge_Peak), pybedtools.BedTool(H3K27ac_Peak), pybedtools.BedTool(chromatin_Peak)


def find_gid(chrom, start, end, PROInfo):
    sidx = int(start)
    eidx = int(end)
    if chrom not in PROInfo:
        return [] 
    overlap_ref = PROInfo[chrom].find(sidx, eidx)
    if not overlap_ref:
        return [] 
    pro_strand = [x.strand for x in overlap_ref]
    if "+" in pro_strand and "-" in pro_strand:
        pro_group = "both"
    elif "+" in pro_strand:
        pro_group = "forward"
    elif "-" in pro_strand:
        pro_group = "reverse"
    else:
        return "."
    return pro_group

def load_peak_yield(mergePeak):
    Peaks = mergePeak.sort().merge(d=0, c='4,5', o="distinct")
    for x in Peaks:
        yield x

def merge_peak_skip_TSS(mergePeak, PromoterInfo, flank = 5000):
    Peaks = load_peak_yield(mergePeak)
    outPeaks = list()
    forward_peak = next(Peaks)
    temp_peak_li = list()
    temp_peak_li.append(list(forward_peak))
    while True:
        try:
            peak_line = next(Peaks)
        except StopIteration:
            break
        chrom, start, end = peak_line[0], int(peak_line[1]), int(peak_line[2])
        forward_chrom, forward_start, forward_end = forward_peak[0], int(forward_peak[1]), int(forward_peak[2])
        distance = start - forward_end
        promoter_info = find_gid(chrom, forward_end, start, PromoterInfo)
        if distance < flank and not promoter_info and forward_chrom==chrom:
            temp_peak_li.append(list(peak_line))
        else:
            temp_merge_peak = pybedtools.BedTool(temp_peak_li).sort().merge(d=flank, c='4,5', o="distinct")
            assert len(temp_merge_peak) == 1
            outPeaks.append(temp_merge_peak[0])
            temp_peak_li = list()
            temp_peak_li.append(list(peak_line)) 
        forward_peak = peak_line
    outPeaks.append(forward_peak)
    return pybedtools.BedTool(outPeaks)

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("--DNase", type=str, dest="DNase", metavar="DNase.peak", required=True)
    base_group.add_argument("--H3K4me4", type=str, dest="H3K4me4", metavar="H3K4me4.peak", required=True)
    base_group.add_argument("--H3K27ac", type=str, dest="H3K27ac", metavar="H3K27ac.peak", required=True)
    base_group.add_argument("--GRO", type=str, dest="GRO", metavar="GRO.peak", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.bed", required=True)
    base_group.add_argument("-prefix", type=str, dest="prefix", metavar="prefix", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_DNase = args.DNase
    f_H3K27ac = args.H3K27ac
    f_H3K4me4 = args.H3K4me4
    f_GRO = args.GRO
    f_enahncer = args.enhancer
    f_promoter = args.promoter
    f_prefix = args.prefix
    f_bed = "/home/wenmiaomiao/cotton/data/Gabor_genome/Gabor.gene.bed6"
    PromoterInfo = write_promoter(f_H3K4me4, f_promoter) 
    PromoterInfo = merge_all_promoter(f_bed, PromoterInfo, flank=2000)
    ChIPPeak = merge_ChIP_peak(f_H3K27ac, PromoterInfo)
    PeakInfo = record_ChIP_DNase(ChIPPeak, f_DNase)

    mergePeak, H3K27acPeak, ChromatinPeak = group_peak(PeakInfo, PromoterInfo)
    H3K27acPeak.saveas(os.path.join(f_prefix,"H3K27ac.Peak"))
    ChromatinPeak.saveas(os.path.join(f_prefix,"Chromatin.Peak"))
    outPeaks = merge_peak_skip_TSS(mergePeak, PromoterInfo, flank=5000)

    GROInfo = load_GRO_info(f_GRO)
    n = 0
    with open(f_enahncer, "w") as fo:
        for x in outPeaks:
            chrom, start, end, group, sample = x[0], x[1], x[2], x[3], x[4]
            start, end = int(start), int(end)
            promoter_info = find_gid(chrom, start, end, PromoterInfo)
            assert not promoter_info
            sample_group = "WT"
            gro_info = find_gid(chrom, start, end, GROInfo)
            if gro_info:
                n = n+1
                peak_name = "enhancer%s" % n    
                outline = [chrom, start, end, peak_name, sample_group, ".", gro_info]
                outline = list(map(str, outline))
                fo.write("\t".join(outline)+"\n")


if __name__ == "__main__":
    main()
