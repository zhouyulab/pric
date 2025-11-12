from collections import defaultdict
import pandas as pd
import sys

def record_total_reads(f_in):
    df = pd.read_csv(f_in, sep="\t")
    df.columns = [x.split("/")[-1].split(".")[0] for x in df.columns]
    read_dict = dict()
    for samp in df.columns[1:]:
        read_dict[samp] = df[samp].sum()
    return read_dict


def main():
    f_RNA = sys.argv[1]
    f_RNA_summary = sys.argv[2]
    f_out = sys.argv[3]
    print(f_out)

    total_read_dict = record_total_reads(f_RNA_summary)
    header = ["Enhancer", "Chrom", "Start", "End", "Strand", "Length","RNA_WT_Rep1", "RNA_WT_Rep2", "RNA_fl_Rep1", "RNA_fl_Rep2", "PRO_WT_Rep1", "PRO_WT_Rep2",]
    with open(f_RNA) as fi, open(f_out, "w") as fo:
        fo.write("\t".join(header)+"\n")
        for line in fi:
            if line.startswith("#"):
                continue
            if line.startswith("Geneid"):
                continue
            data = line.strip().split("\t")
            Enhancer, Chr, Start, End, Strand, Length, WT_Rep1, WT_Rep2, fl_Rep1, fl_Rep2, GRO_Rep1, GRO_Rep2 = data
            WT_Rep1, WT_Rep2, fl_Rep1, fl_Rep2, GRO_Rep1, GRO_Rep2 = int(WT_Rep1), int(WT_Rep2), int(fl_Rep1), int(fl_Rep2), int(GRO_Rep1), int(GRO_Rep2)
            norm_WT_Rep1 = WT_Rep1 * 10000000/total_read_dict["WT_Rep1"]
            norm_WT_Rep2 = WT_Rep2 * 10000000/total_read_dict["WT_Rep2"]
            norm_fl_Rep1 = fl_Rep1 * 10000000/total_read_dict["fl_Rep1"]
            norm_fl_Rep2 = fl_Rep2 * 10000000/total_read_dict["fl_Rep2"]
            norm_GRO_AD_WT_Rep1 = GRO_Rep1 * 10000000/total_read_dict["GRO_AD_WT_Rep1"]
            norm_GRO_AD_WT_Rep2 = GRO_Rep2 * 10000000/total_read_dict["GRO_AD_WT_Rep2"]
            outline = [Enhancer, Chr.split(";")[0], Start.split(";")[1], End.split(";")[1], Strand, Length]
            outline = outline + [norm_WT_Rep1, norm_WT_Rep2, norm_fl_Rep1, norm_fl_Rep2, norm_GRO_AD_WT_Rep1, norm_GRO_AD_WT_Rep2]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()