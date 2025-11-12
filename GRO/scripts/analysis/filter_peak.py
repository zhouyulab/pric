import sys
import pandas as pd
import argparse

def record_total_reads(f_in):
    df = pd.read_csv(f_in, sep="\t")
    df.columns = ["Start", "Rep1", "Rep2"]
    total_reads = [df["Rep1"].sum(), df["Rep2"].sum()]
    return total_reads

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bed", required=True)
    base_group.add_argument("-s", "--summary", type=str, dest="summary", metavar="count.txt.summary", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_in = args.input
    f_out = args.output
    f_summary = args.summary

    total_reads = record_total_reads(f_summary)
    n = 0
    with open(f_in) as fi, open(f_out, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                continue
            if line.startswith("Geneid"):
                continue
            name, chrom, start, end, strand, length, count1, count2= line.strip().split("\t")
            count1, count2 = int(count1), int(count2)
            norm_count1 = count1 * 10000000 / total_reads[0]
            norm_count2 = count2 * 10000000 / total_reads[1]
            if norm_count1 > 10 and norm_count2 > 10:
                n = n + 1
                gro_name = "GRO_%s" % n
                bed_line = [chrom, start, end, gro_name, (norm_count2+norm_count1)/2, strand]
                bed_line = list(map(str, bed_line))
                fo.write("\t".join(bed_line)+"\n")
    

if __name__ == "__main__":
    main()