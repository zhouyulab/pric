import os, sys, argparse
from collections import defaultdict
import pysam
import pandas as pd

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.txt", required=True)
    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fasta", required=True)
    base_group.add_argument("-f", "--flank", type=int, dest="flank", metavar=500, required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_input = args.input
    f_fasta = args.genome
    flank = args.flank
    f_out = args.output

    GWASinfos = pd.read_csv(f_input)
    Fasta = pysam.FastaFile(f_fasta)
    read_indx = 0
    with open(f_out, "w") as f_o:
        for gwas in GWASinfos:
            chrom = gwas.split(":")[-1]
            site = int(gwas.split(":")[-2])
            flank_start, flank_end = max(0, site-flank), site+flank
            seq = Fasta.fetch(chrom, flank_start, flank_end)
            read_indx += 1
            header = "@read%s::%s:%s:%s::%s" % (read_indx,chrom,flank_start,flank_end,site)
            f_o.write(header+"\n")
            f_o.write(seq+"\n")


if __name__ == "__main__":
    main()