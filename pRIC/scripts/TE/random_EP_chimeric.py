import random
import pyBigWig
import sys, argparse
import numpy as np
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

random.seed(234)
def record_chromSize(f_chrom):
    ChromDict = dict()
    with open(f_chrom) as fi:
        for line in fi:
            chrom, size = line.strip().split("\t")
            if "Scaffold" in chrom:
                continue
            ChromDict[chrom] = int(size)
    return ChromDict

def generate_random_site(ChromDict, bw, flank):
    is_overlap = True
    while is_overlap:
        chrom = random.choice(list(ChromDict.keys()))
        site = random.randint(1,ChromDict[chrom])
        start, end = max(site-flank, 1), site + flank
        block = "%s-%s" % (start, end)
        strand = "."
        bw_values = bw.stats(chrom, site-5, site+5, type="mean")
        bw_values = np.nan_to_num(bw_values)
        if not bw_values:
            continue
        block_value = bw_values[0]
        if block_value > 0.1:
            is_overlap = False
            return [chrom, site, strand, block]
        
def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-c", "--chrom", type=str, dest="chrom_size", metavar="chrom.size.txt", required=True)
    base_group.add_argument("-b", "--bw", type=str, dest="bw", metavar="Sample.bw", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("-r", "--round-time", type=int, dest="round_time", metavar=100000, required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_bw = args.bw
    f_chromSize = args.chrom_size
    f_out = args.output
    round_rime = args.round_time

    bw = pyBigWig.open(f_bw)
    ChromDict = record_chromSize(f_chromSize)
    flank = 100
    with open(f_out, "w") as fo:
        for i in range(round_rime):
            donor_chrom, donor_site, donor_strand, donor_block = generate_random_site(ChromDict, bw, flank)
            donor_type = "enhancer"
            donor = "random_enhancer_%s" % (i+1)


            acceptor_chrom, acceptor_site, acceptor_strand, acceptor_block = generate_random_site(ChromDict, bw, flank)
            acceptor_type = "promoter"
            acceptor = "pro_promoter_%s" % (i+1)

            Group = "E-P"

            outline = [donor_chrom, donor_site, donor_strand, donor_block, acceptor_chrom, acceptor_site, acceptor_strand, acceptor_block, 1, donor, acceptor, Group]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()