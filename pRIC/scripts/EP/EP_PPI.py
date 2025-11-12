import os, sys, argparse
from collections import defaultdict
import numpy as np

def record_chrom_size(f_chromsize):
    ChromSize = dict()
    with open(f_chromsize) as f:
        for line in f:
            chrom, size = line.strip().split("\t")
            ChromSize[chrom] = int(size)
    return ChromSize


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="RBP_chimeric_cluster.txt", required=True)
    base_group.add_argument("-e", "--edge", type=str, dest="edge", metavar="edge.tsv", required=True)
    base_group.add_argument("-n", "--node", type=str, dest="node", metavar="note.tsv", required=True)
    base_group.add_argument("--chrom-size", type=str, dest="chrom_size", metavar="chrom.size.", required=True)
    return parser.parse_args(args)

def main():

    args = sys.argv[1:]
    args = parse_args(args)
    f_input = args.input
    f_chromsize = args.chrom_size
    f_node = args.node
    f_edge = args.edge

    ChromSize = record_chrom_size(f_chromsize)
    ChromPair = defaultdict(int)
    with open(f_input) as f_in:
        for line in f_in:
            if line.startswith("Enhancer"):
                continue
            data = line.strip().split("\t")
            enhancer_chrom, promoter_chrom = data[1].split(":")[0], data[3].split(":")[0]
            readName = int(data[4])
            if "Scaffold" in enhancer_chrom or "Scaffold" in promoter_chrom:
                continue
            if enhancer_chrom < promoter_chrom:
                key = (enhancer_chrom, promoter_chrom)
            else:
                key = (promoter_chrom, enhancer_chrom)
            ChromPair[key] += 1

    with open(f_node, "w") as f_n:
        f_n.write("Chrom\tSize\n")
        for chrom in ChromSize:
            if "Scaffold" in chrom:
                continue
            size = ChromSize[chrom]
            f_n.write(chrom + "\t" + str(size) + "\n")

    with open(f_edge, "w") as f_e:
        f_e.write("Node1\tNode2\tNum\tColor\n")
        for chrom1, chrom2 in ChromPair:
            if chrom1.startswith("Ghir_A") and chrom2.startswith("Ghir_A") and chrom1 != chrom2:
                color = "grey70"
            elif chrom1.startswith("Ghir_D") and chrom2.startswith("Ghir_D") and chrom1 != chrom2:
                color = "grey70"
            elif chrom1 == "Ghir_A01":
                color = "A01:#DF2135"
            elif chrom1 == "Ghir_A02":
                color = "A02:#FF33CC"
            elif chrom1 == "Ghir_A03":
                color = "#A03:F6BA4A"
            elif chrom1 == "Ghir_A04":
                color = "A04:#FFA28B"
            elif chrom1 == "Ghir_A05":
                color = "A05:#FC6882"
            elif chrom1 == "Ghir_A06":
                color = "A06:DF65B0"
            elif chrom1 == "Ghir_A07":
                color = "A07:#C994C7"
            elif chrom1 == "Ghir_A08":
                color = "A08:#008856"
            elif chrom1 == "Ghir_A09":
                color = "A09:#41AB5D"
            elif chrom1 == "Ghir_A10":
                color = "A10:#007BC3"
            elif chrom1 == "Ghir_A11":
                color = "A11:#4292C6"
            elif chrom1 == "Ghir_A12":
                color = "A12:#7EB2E4"
            elif chrom1 == "Ghir_A13":
                color = "A13:#FFCC33"

            num = ChromPair[(chrom1, chrom2)]
            outline = [chrom1, chrom2, str(num), color]
            f_e.write("\t".join(outline)+"\n")


if __name__ == "__main__":
    main()
