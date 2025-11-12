from scipy.stats import fisher_exact, false_discovery_control
from collections import defaultdict
import sys,os
import argparse

def record_EP_reads(f_in):
    total_reads = 0
    featureCount = defaultdict(int)
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split("\t")
            feature1, feature2, reads = data[3], data[9], int(data[-1])
            total_reads += reads
            featureCount[feature1] += reads
            featureCount[feature2] += reads
    return featureCount, total_reads

def calculate_pvalue(f_in, featureCount, total_reads):
    pvalue_li = []
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split("\t")
            feature1, feature2, read_num = data[3], data[9], int(data[-1])
            feature1_num = featureCount[feature1]
            feature2_num = featureCount[feature2]
            odds, pvalue = fisher_exact([[read_num, feature1_num-read_num], [feature2_num-read_num, total_reads - feature1_num-feature2_num+read_num]])
            pvalue_li.append(pvalue)
    FDR_li = false_discovery_control(pvalue_li, method='bh')
    return pvalue_li, FDR_li


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="EP.input", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="EP.output", required=True)
    base_group.add_argument("--out-EP", type=str, dest="output_EP", metavar="EP.output", required=True)
    base_group.add_argument("-c", "--cutoff", type=int, dest="cutoff", metavar="EP.output", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_input = args.input
    f_output = args.output
    f_outEP = args.output_EP
    reads_cutoff = args.cutoff

    header = ["Feature1", "Feature1Pos", "Feature2", "Feature2Pos", "Reads", "Group", "p-value", "FDR"]
    header_EP = ["Enhancer", "EnhancerPos", "Promoter", "PromoterPos", "Reads", "Group", "p-value", "FDR"]

    featureCount, total_reads = record_EP_reads(f_input)
    pvalue_li, FDR_li = calculate_pvalue(f_input, featureCount, total_reads)
    with open(f_input) as fi, open(f_output, "w") as fo, open(f_outEP, "w") as f_ep:
        fo.write("\t".join(header)+"\n")
        f_ep.write("\t".join(header_EP)+"\n")
        for indx,line in enumerate(fi):
            data = line.strip().split("\t")
            feature1, feature2, read_num = data[3], data[9], int(data[-1])
            feature1_pos = "%s:%s-%s" % (data[0], data[1], data[2])
            feature2_pos = "%s:%s-%s" % (data[6], data[7], data[8])
            pvalue = pvalue_li[indx]
            FDR = FDR_li[indx]
            if feature1 == feature2:
                continue
            if read_num < reads_cutoff:
                continue
            if pvalue <0.05 and FDR < 0.05:
                if feature1.startswith("pro") and feature2.startswith("pro"):
                    interaction_group = "P-P"
                    outline = [feature1, feature1_pos, feature2, feature2_pos, read_num, interaction_group, pvalue, FDR]
                elif feature1.startswith("pro") and not feature2.startswith("pro"):
                    interaction_group = "E-P"
                    outline = [feature2, feature2_pos, feature1, feature1_pos, read_num, interaction_group, pvalue, FDR]
                    outline = list(map(str, outline))
                    f_ep.write("\t".join(outline)+"\n")
                elif not feature1.startswith("pro") and feature2.startswith("pro"):
                    interaction_group = "E-P"
                    outline = [feature1, feature1_pos, feature2, feature2_pos, read_num, interaction_group, pvalue, FDR]
                    outline = list(map(str, outline))
                    f_ep.write("\t".join(outline)+"\n")
                else:
                    interaction_group = "E-E"    
                    outline = [feature1, feature1_pos, feature2, feature2_pos, read_num, interaction_group, pvalue, FDR]
                outline = list(map(str, outline))
                fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()

