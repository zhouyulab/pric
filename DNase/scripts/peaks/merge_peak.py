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

    base_group.add_argument("--input-rep1", type=str, dest="input_rep1", metavar="DNase_rep1.narrowPeak", required=True)
    base_group.add_argument("--input-rep2", type=str, dest="input_rep2", metavar="DNase_rep2.narrowPeak", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="outdir", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_rep1 = args.input_rep1
    f_rep2 = args.input_rep2
    f_out = args.output

    Peak_Rep1 = pybedtools.BedTool(f_rep1).merge(d=200)
    Peak_Rep2 = pybedtools.BedTool(f_rep2).merge(d=200)

    fo = open(f_out, "w")
    merge_Peak = Peak_Rep1.intersect(Peak_Rep2, wo=True)
    n = 0
    for x in merge_Peak:
        merge_start = min(x[1],x[4])
        merge_end = max(x[2], x[5])
        n = n + 1
        name = "WT_DNase_%s" % (n)
        outline = [x[0], merge_start, merge_end, "WT", name, "."]
        outline = list(map(str, outline))
        fo.write("\t".join(outline)+"\n")

    fo.close()

if __name__ == "__main__":
    main()