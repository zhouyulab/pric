import pybedtools
import sys
import argparse


def outGTF(chrom, start, end, strand, proName, handle):
        gene_info = 'gene_type "%s"; gene_name "%s";'%("GRO", proName)
        trans_info = 'transcript_id "%s"; gene_type "%s"; gene_name "%s";'%(proName, "GRO", proName)
        gtf_info_1 = [chrom, "GRO", "gene", start, end, ".", strand, ".", gene_info]
        gtf_info_2 = [chrom, "GRO", "transcript", start, end, ".", strand, ".", trans_info]
        gtf_info_3 = [chrom, "GRO", "exon", start, end, ".", strand, ".", trans_info]
        for info in [gtf_info_1, gtf_info_2, gtf_info_3]:
            handle.write("\t".join(info) + "\n")

def merge_rep(f_rep1, f_rep2):
    Peaks = list()
    with open(f_rep1) as fi:
        for line in fi:
            data = line.strip().split("\t")
            Peaks.append(data)
    with open(f_rep2) as fi:
        for line in fi:
            data = line.strip().split("\t")
            Peaks.append(data)
    Peaks = pybedtools.BedTool(Peaks)
    Peaks_merge = Peaks.sort().merge(d=0)
    return Peaks_merge

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("--prak-R-1", type=str, dest="peak_Reverse_rep1", metavar="Forward_peaks.narrowPeak", required=True)
    base_group.add_argument("--peak-F-1", type=str, dest="peak_Forward_rep1", metavar="Reverse_peaks.narrowPeak", required=True)
    base_group.add_argument("--prak-R-2", type=str, dest="peak_Reverse_rep2", metavar="Forward_peaks.narrowPeak", required=True)
    base_group.add_argument("--peak-F-2", type=str, dest="peak_Forward_rep2", metavar="Reverse_peaks.narrowPeak", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed", required=True)
    base_group.add_argument("-gtf", "--gtf", type=str, dest="GTF", metavar="output.gtf", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_forward_rep1 = args.peak_Forward_rep1
    f_reverse_rep1 = args.peak_Reverse_rep1
    f_forward_rep2 = args.peak_Forward_rep2
    f_reverse_rep2 = args.peak_Reverse_rep2
    f_bed = args.output
    f_gtf = args.GTF


    peak_Forward = merge_rep(f_forward_rep1, f_forward_rep2)
    peak_Reverse = merge_rep(f_reverse_rep1, f_reverse_rep2)

    with open(f_bed, "w") as fo, open(f_gtf, "w") as handle:
        n = 0
        for x in peak_Forward:
            chrom = x[0]
            start = str(int(x[1])+1)
            end = str(int(x[2])+1)
            strand = "-"
            n = n+1
            name = "GRO%s" % n
            outline = [chrom, str(start), str(end), name, ".", strand]
            fo.write("\t".join(outline)+"\n")
            outGTF(chrom, start, end, strand, name, handle)

        for y in peak_Reverse:
            chrom = y[0]
            start = str(int(y[1])+1)
            end = str(int(y[2])+1)
            strand = "+"
            n = n+1
            name = "GRO%s" % n
            outline = [chrom, str(start), str(end), name, ".", strand]
            fo.write("\t".join(outline)+"\n")
            outGTF(chrom, start, end, strand, name, handle)

if __name__ == "__main__":
    main()