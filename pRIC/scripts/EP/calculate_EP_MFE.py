import random
import pysam
import os, sys, argparse
from Bio.Seq import Seq
from collections import defaultdict
random.seed(234)

def calculate_MFE(enhancer_seq, promoter_seq, enhancer_key, promoter_key, f_outdir, prefix):
    f_path = os.path.join(f_outdir, enhancer_key+"__"+promoter_key)
    if not os.path.exists(f_path):
        os.mkdir(f_path)

    f_enhancer = os.path.join(f_outdir, enhancer_key+"__"+promoter_key, "%s.enhancer.fasta" % prefix)
    with open(f_enhancer, "w") as f_e:
        f_e.write(">"+enhancer_key +"\n")
        f_e.write(enhancer_seq)

    f_promoter = os.path.join(f_outdir, enhancer_key+"__"+promoter_key,  "%s.promoter.fasta" % prefix)
    with open(f_promoter, "w") as f_p:
        f_p.write(">"+promoter_key+"\n")
        f_p.write(promoter_seq)
    
    f_cnt = os.path.join(f_outdir, enhancer_key+"__"+promoter_key, "%s.cnt.txt" % prefix)
    command1 = "DuplexFold -d %s %s %s" % (f_enhancer, f_promoter, f_cnt)
    os.system(command1)

    f_ct = os.path.join(f_outdir, enhancer_key+"__"+promoter_key, "%s.ct" % prefix)
    command2 = "ct2dot %s 1 %s" % (f_cnt, f_ct)
    os.system(command2)

    with open(f_cnt) as f:
        for line in f:
            if "ENERGY" in line:
                MFE = line.strip().split("  ")[1].split("=")[1].strip()
    return MFE

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fa", required=True)
    base_group.add_argument("-t", "--tem_dir", type=str, dest="tem_dir", metavar="promoter.txt", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="chimeric", metavar="chimeric.txt", required=True)
    base_group.add_argument("--flank", type=int, dest="flank", metavar=200, required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_fasta = args.genome
    f_chimeric = args.chimeric
    f_outdir = args.tem_dir
    f_out = args.output
    flank = args.flank

    header = ["Enhancer", "EnhancerChrom", "EnhancerSite", "Promoter", "PromoterChrom", "PromoterSite", "Group", "MFE"]
    Fasta = pysam.FastaFile(f_fasta)

    storage_pair = set()
    with open(f_chimeric) as fi, open(f_out, "w") as fo:
        fo.write("\t".join(header)+"\n")
        for line in fi:
            if line.startswith("DonorChrom"):
                continue
            data = line.strip().split("\t")
            DonorChrom, DonorCJS, DonorStrand, DonorBlocks, AcceptorChrom, AcceptorCJS, AcceptorStrand, AcceptorBlocks, ReadNum,DonorFeature,AcceptorFeature, Group =data
            if Group != "E-P":
                continue
            if DonorFeature.startswith("pro"):
                assert not AcceptorFeature.startswith("pro")
                promoter_chrom, promoter_site, promoter_strand, promoter = DonorChrom, int(DonorCJS), DonorStrand, DonorFeature
                enhancer_chrom, enhancer_site, enhancer_strand, enhancer = AcceptorChrom, int(AcceptorCJS), AcceptorStrand, AcceptorFeature

            else:
                assert AcceptorFeature.startswith("pro")
                assert not DonorFeature.startswith("pro")
                promoter_chrom, promoter_site, promoter_strand, promoter = AcceptorChrom, int(AcceptorCJS), AcceptorStrand, AcceptorFeature
                enhancer_chrom, enhancer_site, enhancer_strand, enhancer = DonorChrom, int(DonorCJS), DonorStrand, DonorFeature
            
            key = (enhancer, enhancer_chrom, enhancer_site, promoter, promoter_chrom, promoter_site)
            if key in storage_pair:
                continue

            storage_pair.add(key)
            promoter_start, promoter_end = max(1, promoter_site-flank), promoter_site+flank
            promoter_seq = Fasta.fetch(promoter_chrom, promoter_start, promoter_end)
            promoter_seq = Seq(promoter_seq)
            if promoter_strand=="-":
                promoter_seq = promoter_seq.reverse_complement()
            promoter_rna = promoter_seq.transcribe()
            promoter_key = "%s_%s" % (promoter, promoter_site)

            enhancer_start, enhancer_end = max(1, enhancer_site-flank), enhancer_site+flank
            enhancer_seq = Fasta.fetch(enhancer_chrom, enhancer_start, enhancer_end)
            enhancer_seq = Seq(enhancer_seq)
            if enhancer_strand=="-":
                enhancer_seq = enhancer_seq.reverse_complement()
            enhancer_rna = enhancer_seq.transcribe()
            enhancer_key = "%s_%s" % (enhancer, enhancer_site)

            EP_MFE = calculate_MFE(str(enhancer_rna), str(promoter_rna), enhancer_key, promoter_key, f_outdir, "EP")

            promoter_shuffle_rna = list(str(promoter_rna))
            random.shuffle(promoter_shuffle_rna)
            enhancer_shuffle_rna = list(str(enhancer_rna))
            random.shuffle(enhancer_shuffle_rna)
            shuffle_MFE = calculate_MFE("".join(enhancer_shuffle_rna), "".join(promoter_shuffle_rna), enhancer_key, promoter_key, f_outdir, "shuffle")

            outline1 = [enhancer, enhancer_chrom, enhancer_site, promoter, promoter_chrom, promoter_site, "EP", EP_MFE]
            outline1 = list(map(str, outline1))
            fo.write("\t".join(outline1)+"\n")

            outline2 = [enhancer, enhancer_chrom, enhancer_site, promoter, promoter_chrom, promoter_site, "shuffle", shuffle_MFE]
            outline2 = list(map(str, outline2))
            fo.write("\t".join(outline2)+"\n")

if __name__ == "__main__":
    main()
