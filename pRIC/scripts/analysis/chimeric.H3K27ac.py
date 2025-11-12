import argparse
import os
import sys
from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import numpy as np

class PeakLine(object):
    def __init__(self, line):
        # print(line)
        data = line.strip().split("\t")
        if len(data) == 10:
            chrom, start, end, name, score, strand, signalvalue, pvalue, pvalue, peak_center = data
            self.peak_center = int(start) + int(peak_center)
        elif len(data) == 6:
            chrom, start, end, sample, name, strand = data
            self.peak_center = (int(start) + int(end)) /2
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name


    def get_indx(self, site):
        if site > self.peak_center:
            indx = 1000 + (site-self.peak_center)
        else:
            indx = 1000 - (self.peak_center-site)
        return indx

def load_narrowPeak(f_peak):
    PeakInfo = dict()
    with open(f_peak) as f:
        for line in f:
            peak_line = PeakLine(line)
            data = line.strip().split("\t")
            chrom, start, end = data[:3]
            key = chrom
            if key not in PeakInfo:
                PeakInfo[key] = IntervalTree()
            PeakInfo[key].insert(int(start), int(end), peak_line)
    return PeakInfo


def main():
    parser = argparse.ArgumentParser(description="Fetch gene interactions")
    parser.add_argument("-peak", dest="peak", help="H3K27ac_peaks.narrowPeak")
    parser.add_argument("-c", dest="chimeric", help="WT_ChimericGroup.txt")
    parser.add_argument("-o", dest="output", help="output.tsv")
    args = parser.parse_args()
    print(args)
    f_peak = args.peak
    f_chimeric = args.chimeric
    f_out = args.output

    PeakInfo = load_narrowPeak(f_peak)
    PeakValue = dict()
    with open(f_chimeric) as fi:
        for line in fi:
            data = line.strip().split("\t")
            if data[0] == "DonorChrom":
                continue
            DonorChrom,DonorCJS,DonorStrand,DonorBlocks,AcceptorChrom,AcceptorCJS,AcceptorStrand,AcceptorBlocks,ReadNum,DonorGene,AccetperGene,Group = data
            ReadNum = int(ReadNum)
            DonorCJS, AcceptorCJS = int(DonorCJS), int(AcceptorCJS)
            if DonorGene == "None" and AccetperGene != "None":
                if DonorChrom in PeakInfo.keys():
                    overlap_peak = PeakInfo[DonorChrom].find(DonorCJS, DonorCJS+1)
                    for peak in overlap_peak:
                        if peak.name not in PeakValue:
                            PeakValue[peak.name] = np.zeros(2000)
                        chimeric_indx = peak.get_indx(DonorCJS)
                        chimeric_indx = int(chimeric_indx)
                        if chimeric_indx > 0 and chimeric_indx < 2000:
                            PeakValue[peak.name][chimeric_indx] += ReadNum

            if AccetperGene == "None" and DonorGene!= "None":
                if AcceptorChrom in PeakInfo.keys():
                    overlap_peak = PeakInfo[AcceptorChrom].find(AcceptorCJS, AcceptorCJS+1)
                    for peak in overlap_peak:
                        if peak.name not in PeakValue:
                            PeakValue[peak.name] = np.zeros(2000)
                        chimeric_indx = peak.get_indx(AcceptorCJS)
                        chimeric_indx = int(chimeric_indx)
                        if chimeric_indx > 0 and chimeric_indx < 2000:
                            PeakValue[peak.name][chimeric_indx] += ReadNum

    with open(f_out, "w") as fo:
        fo.write("Index"+"\t"+"Reads"+"\t"+ "Peak"+"\n")
        for i in range(2000):
            indx_value = 0
            for peak_name in PeakValue:
                indx_value += PeakValue[peak_name][i]
            fo.write(str(i)+"\t"+str(indx_value)+"\t"+"H3K27ac"+"\n")

if __name__ == "__main__":
    main()


