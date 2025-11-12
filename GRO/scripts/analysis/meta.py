import argparse as ap
from os import SCHED_OTHER
import numpy as np
import pyBigWig

class Bed3(object):
    def __init__(self, chrom, chromStart, chromEnd, strand):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.strand = strand
    
    def __str__(self):
        return "%s:%d-%d" % self.change_to_list()[:3]
    
    def change_to_list(self):
        """
        Get a list of range.
        """
        return tuple(self.__dict__.values())

    def fetch_bw(self, bw):
        if self.chromStart == self.chromEnd:
            return 0
        else:
            bw_values = bw.values(*self.change_to_list())
            cov = np.mean(np.nan_to_num(bw_values))
            return np.abs(cov)

class Bed6(object):
    def __init__(self, record_line):
        record_list = record_line.strip().split('\t')
        self.chrom = record_list[0]
        self.chromStart = int(record_list[1])
        self.chromEnd = int(record_list[2])
        self.name = record_list[3]
        self.score = record_list[4]
        self.strand = record_list[5]

    def give_upstream(self, initial, length, chromDict):
        if self.strand == '+':
            site = max(1, initial - length)
        elif self.strand == '-':
            site = min(initial + length + 1, chromDict[self.chrom])
        return site

    def give_downstream(self, initial, length, chromDict):
        if self.strand == '+':
            site =  min(initial + length + 1, chromDict[self.chrom])
        elif self.strand == '-':
            site = max(1, initial - length)
        return site

    def choose_bw(self, bw1, bw2):
        """
        If the extand range is positive strand, use bw1, else use bw2.
        """
        bw = bw1 if self.strand == '+' else bw2
        return bw

    def slide_windows(self, before, after, bin_num_up, bin_num_body, bin_num_down, chromDict):
        """
        From upstream to downstream.
        Distinguish forward and reverse strand.
        Scale region model.
        Calculate extended region and scaled body region respectively.
        """
        if self.strand == "+":
            RegionStart = self.give_upstream(self.chromStart, before, chromDict)
            RegionEnd = self.give_downstream(self.chromEnd, after, chromDict)
            upstream_list = np.linspace(RegionStart, self.chromStart, bin_num_up + 1, endpoint=True, retstep=False, dtype=int)
            downstream_list = np.linspace(self.chromEnd, RegionEnd, bin_num_down+ 1, endpoint=True, retstep=False, dtype=int)
        if self.strand == '-':
            RegionStart = self.give_downstream(self.chromStart, after, chromDict)
            RegionEnd = self.give_upstream(self.chromEnd, before, chromDict)
            upstream_list = np.linspace(self.chromEnd, RegionEnd, bin_num_up + 1, endpoint=True, retstep=False, dtype=int)
            downstream_list = np.linspace(RegionStart, self.chromStart, bin_num_down+ 1, endpoint=True,retstep=False,dtype=int)

        body_list = np.linspace(self.chromStart, self.chromEnd, bin_num_body + 1, endpoint=True,retstep=False,dtype=int)

        up_slide_list = []
        body_slide_list = []
        down_slide_list = []

        if self.strand == '+':
            for i in range(bin_num_up):
                i_start = upstream_list[i]
                i_end = upstream_list[i+1]
                up_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))
            for i in range(bin_num_body):
                i_start = body_list[i]
                i_end = body_list[i+1]
                body_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))
            for i in range(bin_num_down):
                i_start = downstream_list[i]
                i_end = downstream_list[i+1]
                down_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))               
        if self.strand == '-':
            for i in range(bin_num_up, 0, -1):
                i_start = upstream_list[i-1]
                i_end = upstream_list[i]
                up_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))
            for i in range(bin_num_body, 0, -1):
                i_start = body_list[i-1]
                i_end = body_list[i]
                body_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))
            for i in range(bin_num_down, 0, -1):
                i_start = downstream_list[i-1]
                i_end = downstream_list[i]
                down_slide_list.append(Bed3(self.chrom, i_start, i_end, self.strand))

        slide_list = up_slide_list + body_slide_list + down_slide_list + [self.name]
        return slide_list

    def output_key(self):
        return [self.chrom, self.chromStart, self.chromEnd, self.name, self.strand]

def record_chrom(f_chrom_size):
    chrom_size_dict = dict()
    with open(f_chrom_size) as fi:
        for line in fi:
            data = line.strip().split("\t")
            chrom_size_dict[data[0]] = int(data[1])
    return chrom_size_dict

def get_slide_value(slide_list, bwR, bwF):
    template_slide_cov = []
    antisense_slide_cov  = []
    for i, s1, in enumerate(slide_list):
        if s1.strand == "+":
            template_slide_cov.append(s1.fetch_bw(bwF))
            antisense_slide_cov.append(s1.fetch_bw(bwR))
        else:
            template_slide_cov.append(s1.fetch_bw(bwR))
            antisense_slide_cov.append(s1.fetch_bw(bwF))
    return template_slide_cov, antisense_slide_cov

if __name__ == "__main__":
    p = ap.ArgumentParser(
        description = __doc__,
        formatter_class = ap.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--input", type=str, required=True, help="bed file")
    p.add_argument("--bw-R", type=str, required=True, help="neg strand bigwig")
    p.add_argument("--bw-F", type=str, required=True, help="pos strand bigwig")
    p.add_argument("--output", type=str, required=True, help="output coverage")
    p.add_argument("--after", type=int, required=True, help="downstream width")
    p.add_argument("--before", type=int, required=True, help="upstream width")
    p.add_argument("--around", type=int, required=True, help="gene body exclued width")
    p.add_argument("--bin_num_up", type=int, required=True, help="bin number")
    p.add_argument("--bin_num_body", type=int, required=True, help="bin number")
    p.add_argument("--bin_num_down", type=int, required=True, help="bin number")
    p.add_argument("--chrom_size", type=str, required=True, help="chrom size")
    args = p.parse_args()

    chromSize_dict = record_chrom(args.chrom_size)
    with open(args.input, 'r') as fi, open(args.output, 'w') as fo:
        bwR = pyBigWig.open(args.bw_R)
        bwF = pyBigWig.open(args.bw_F)
        for record_line in fi:
            input_region = Bed6(record_line)
            if input_region.chrom.startswith("Scaffold") or (input_region.chromEnd - input_region.chromStart) <= args.bin_num_body:
                continue

            slide_all = input_region.slide_windows(before=args.before, after=args.after,
                         bin_num_up=args.bin_num_up,bin_num_body=args.bin_num_body, bin_num_down=args.bin_num_down, 
                         chromDict = chromSize_dict)
            slide_list = slide_all[:-1]
            gene = slide_all[-1]
            

            template_slide_cov, antisense_slide_cov = get_slide_value(slide_list, bwR, bwF)
                
         
            bins_info_template = [','.join(map(str, template_slide_cov))]
            bins_info_antisense = [','.join(map(str, antisense_slide_cov))]
            outline_template = [gene, "template"] + bins_info_template
            outline_antisense = [gene, "antisense"] + bins_info_antisense
            print(*outline_template, sep=',', file=fo)
            print(*outline_antisense, sep=',', file=fo)
