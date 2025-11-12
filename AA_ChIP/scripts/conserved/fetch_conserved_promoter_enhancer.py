import os, sys, argparse
from collections import defaultdict
import pybedtools

def filter_blast(f_in, f_out):
    with open(f_in) as fi, open(f_out, "w") as fo:
        for line in fi:
            tetraploid, diploid, identity, alignment_length, mismatches, gaps, tetraploid_start, tetraploid_end, diploid_start, diploid_end,  e_value, score = line.strip().split("\t")
            alignment_length, identity = int(alignment_length), float(identity)

            diploid_feature, diploid_pos = diploid.split("::")
            diploid_start, diploid_end = diploid_pos.split(":")[1].split("-")
            diploid_length = abs(int(diploid_end)-int(diploid_start))
            diploid_ratio = alignment_length/diploid_length

            tetraploid_feature, tetraploid_pos = tetraploid.split("::")
            tetraploid_start, tetraploid_end = tetraploid_pos.split(":")[1].split("-")
            tetraploid_length = abs(int(tetraploid_end)-int(tetraploid_start))
            tetraploid_ratio = alignment_length/tetraploid_length

            # if tetraploid_feature == "E8754":
            #     print(tetraploid_feature, diploid_feature, identity, diploid_ratio, tetraploid_ratio)
            if identity > 90 and max(diploid_ratio, tetraploid_ratio) > 0.6:
                outline = [diploid_feature, diploid_pos, tetraploid_feature, tetraploid_pos, alignment_length, identity, diploid_ratio, tetraploid_ratio]
                outline = list(map(str, outline))
                fo.write("\t".join(outline)+"\n")


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.bed", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("--out-promoter", type=str, dest="output_promoter", metavar="promoter.pair.txt", required=True)
    base_group.add_argument("--out-enhancer", type=str, dest="output_enhancer", metavar="enhancer.pair.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_promoter = args.promoter
    f_enhancer = args.enhancer
    f_output_enhancer = args.output_enhancer
    f_output_promoter = args.output_promoter

    filter_blast(f_promoter, f_output_promoter)
    filter_blast(f_enhancer, f_output_enhancer)

if __name__ == "__main__":
    main()