import argparse
import sys
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="ChimericRead.tsv", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    return parser.parse_args(args)


def cnt_reads(f_in, f_out):
    cnt_dict = defaultdict(int)
    total_reads = 0
    with open(f_in, "r") as f:
        for indx, line in enumerate(f.readlines()):
            if indx == 0:
                continue
            data = line.rstrip("\n").split("\t")
            total_reads = total_reads + 1
            DonorChrom, DonorJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorJS, AcceptorStrand, AcceptorBlocks, ReadName = data
            key = (DonorChrom, DonorJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorJS, AcceptorStrand, AcceptorBlocks)
            cnt_dict[key] += 1

    with open(f_out, "w") as f:
        header = "\t".join(["DonorChrom", "DonorCJS", "DonorStrand", "DonorBlocks", "AcceptorChrom", "AcceptorCJS", "AcceptorStrand", "AcceptorBlocks", "ReadNum"]) + "\n"
        f.write(header)
        for key, cnt in sorted(cnt_dict.items()):
            data = list(key) + [str(cnt)]
            f.write("\t".join(data)+"\n")


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    cnt_reads(f_in, f_out)


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
