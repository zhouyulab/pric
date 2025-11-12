import os, sys, argparse
from collections import defaultdict
import pybedtools

def calculate_distance(enhancer, promoter):
    enahncer_start = int(enhancer[1])
    enahncer_end = int(enhancer[2])
    promoter_start = int(promoter[1])
    promoter_end = int(promoter[2])
    if promoter_start >= enahncer_end:
        distance = promoter_start - enahncer_end
        return ["right", distance]
    else:
        distance = enahncer_start - promoter_end
        return ["left", distance]

def record_promoter(f_promoter):
    Promoter = dict()
    with open(f_promoter) as f:
        for line in f:
            data = line.strip().split("\t")
            if data[0] not in Promoter:
                Promoter[data[0]] = dict()
            Promoter[data[0]][data[3]] = data
    return Promoter



def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.bed", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="outdir", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_promoter = args.promoter
    f_enhancer = args.enhancer
    f_output = args.output

    Promoter = record_promoter(f_promoter)
    with open(f_enhancer) as f, open(f_output, "w") as f_o:
        for line in f:
            data = line.strip().split("\t")
            if data[0] not in Promoter:
                outline = data + ["None", "Inf"]
                f_o.write("\t".join(outline)+"\n")
                continue
            promoter_li = Promoter[data[0]]
            left_distance_li = list()
            right_distance_li = list()
            for p in promoter_li:
                pos, diatance = calculate_distance(data, promoter_li[p])
                if pos == "left":
                    left_distance_li.append([p, diatance])
                if pos == "right":
                    right_distance_li.append([p, diatance])
        
            if left_distance_li:
                left_nearst_promoter = sorted(left_distance_li, key=lambda x: x[1])[0][0]
            else:
                left_nearst_promoter = "NA"
            if right_distance_li:
                right_nearst_promoter = sorted(right_distance_li, key=lambda x: x[1])[0][0]
            else:
                right_nearst_promoter = "NA"
            outline = data + [left_nearst_promoter, right_nearst_promoter]
            f_o.write("\t".join(outline)+"\n")
# 

if __name__ == "__main__":
    main()