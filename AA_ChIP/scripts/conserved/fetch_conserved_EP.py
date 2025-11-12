from collections import defaultdict
import os, sys, argparse

def record_pair(f_in):
    PairInfo = defaultdict(list)
    with open(f_in) as fi:
        for line in fi:
            data = line.strip().split("\t")
            diploid, tetraploid = data[0], data[2]
            PairInfo[tetraploid].append(diploid)
    return PairInfo

def record_EP(f_in):
    PairInfo = defaultdict(list)
    with open(f_in) as fi:
        for line in fi:
            if line.startswith("Enhancer"):
                continue
            data = line.strip().split("\t")
            enhancer, promoter = data[0], data[2]
            PairInfo[enhancer].append(promoter)
    return PairInfo

def get_conserved_promoter(PromoterPair, promoter_li):
    conserved_promoter = list()
    for promoter in promoter_li:
        temp_promoter = PromoterPair[promoter]
        conserved_promoter.extend(temp_promoter)
    conserved_promoter = list(set(conserved_promoter))
    return conserved_promoter


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--promoter", type=str, dest="promoter", metavar="promoter.bed", required=True)
    base_group.add_argument("-e", "--enhancer", type=str, dest="enhancer", metavar="enhancer.bed", required=True)
    base_group.add_argument("--diploid", type=str, dest="diploid_EP", metavar="A_WT_merge.EP.fisher.txt", required=True)
    base_group.add_argument("--tetraploid", type=str, dest="tetraploid_EP", metavar="AD_pCp_merge.EP.fisher.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_promoter_pair = args.promoter
    f_enahcner_pair = args.enhancer
    f_tetraploid_EP = args.tetraploid_EP
    f_diploid_EP = args.diploid_EP
    f_out = args.output


    PromoterPair = record_pair(f_promoter_pair)
    EnhancerPair = record_pair(f_enahcner_pair)

    TetraploidEP = record_EP(f_tetraploid_EP)
    DiploidEP = record_EP(f_diploid_EP)


    header = ["AD.enhancer", "AD.num", "AD.Promoter", "AD.PromoterPair", "A.enhancer", "A.num", "A.Promoter", "Group"]
    with open(f_out, "w") as fo:
        fo.write("\t".join(header)+"\n")
        for tetra_enhancer in TetraploidEP:
            di_enahncer_li = EnhancerPair[tetra_enhancer]
            if not di_enahncer_li:
                continue
            for di_enahncer in di_enahncer_li:
                tetra_promoter = TetraploidEP[tetra_enhancer]
                tetra_promoter = list(set(tetra_promoter))
                di_promoter = DiploidEP[di_enahncer]
                if not di_promoter:
                    continue
                tetra_promoter_pair = get_conserved_promoter(PromoterPair, tetra_promoter)
                group = "conserved" if (set(di_promoter) & set(tetra_promoter_pair)) else "novel"
                outline = [tetra_enhancer, len(tetra_promoter), ";".join(tetra_promoter), ";".join(tetra_promoter_pair)]
                outline = outline + [di_enahncer, len(di_promoter), ";".join(di_promoter), group]

                outline = list(map(str, outline))
                fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()

# python scripts/conserved/fetch_conserved_EP.py -p ../results/conserved/promoter.pair.txt \
# -e ../results/conserved/enhancer.pair.txt \
# --diploid ../../AA_RIC-seq/results2/EP2/A_WT_merge/A_WT_merge.EP.fisher.txt \
# --tetraploid ../../RIC-seq/results2/EP2/AD_pCp_merge/AD_pCp_merge.EP.fisher.txt \
# -o ../results/conserved/EP.conserved.pair.txt