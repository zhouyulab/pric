import pybedtools
import sys
import argparse

def load_peak(f_in, f_out):
    Peaks = list()
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split()
            # print(data)
            chrom, start, end, name, num1, strand, fold_enrichment, log10_pvalue, log10_qvalue, num2 = data
            name_prefix = "_".join(name.split("_")[:-1])
            start, end = int(start), int(end)
            fold_enrichment = float(fold_enrichment)
            if fold_enrichment <= 5:
                continue
            peak_line = [chrom, start, end, str(fold_enrichment), log10_pvalue, strand]
            Peaks.append(peak_line)
    Peaks = pybedtools.BedTool(Peaks)
    Peaks = Peaks.merge(d=200)

    n = 0
    peakInfo = dict()
    with open(f_out, "w") as fo:
        for x in Peaks:
            n += 1
            peak_name = "%s%s" % (name_prefix, n)
            peakInfo[tuple(x)] = peak_name 
            outline = list(x) + [peak_name, ".", "."]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")
    return Peaks, peakInfo


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("--input-WT", type=str, dest="input_WT", metavar="WT_H3K27ac_peaks.narrowPeak", required=True)
    base_group.add_argument("--input-fl", type=str, dest="input_fl", metavar="fl_H3K27ac_peaks.narrowPeak", required=True)
    base_group.add_argument("--output-WT", type=str, dest="output_WT", metavar="WT_H3K27ac.merge.Peak", required=True)
    base_group.add_argument("--output-fl", type=str, dest="output_fl", metavar="fl_H3K27ac.merge.Peak", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="outdir", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_WT_input = args.input_WT
    f_WT_output = args.output_WT
    f_fl_input = args.input_fl
    f_fl_output = args.output_fl
    f_out = args.output

    WT_Peak, WTInfo = load_peak(f_WT_input, f_WT_output)
    fl_Peak, flInfo = load_peak(f_fl_input, f_fl_output)

    fo = open(f_out, "w")
    name_set = set()
    merge_Peak = WT_Peak.intersect(fl_Peak, wo=True)
    for x in merge_Peak:
        wt_key = (x[0], x[1], x[2])
        fl_key = (x[3], x[4], x[5])
        wt_name = WTInfo[wt_key]
        fl_name = flInfo[fl_key]
        name_set.add(wt_name)
        name_set.add(fl_name)
        merge_start = min(x[1],x[4])
        merge_end = max(x[2], x[5])
        outline = [x[0], merge_start, merge_end, "WT,fl", "%s,%s" % (wt_name, fl_name), "."]
        outline = list(map(str, outline))
        fo.write("\t".join(outline)+"\n")

    for y in WT_Peak:
        wt_name = WTInfo[tuple(y)]
        if wt_name not in name_set:
            outline = list(y) + ["WT", wt_name, "."]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

    for y in fl_Peak:
        fl_name = flInfo[tuple(y)]
        if fl_name not in name_set:
            outline = list(y) + ["fl", fl_name, "."]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

    fo.close()

if __name__ == "__main__":
    main()