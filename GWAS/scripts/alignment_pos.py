import re, sys
import pysam
import argparse
from collections import defaultdict

def record_refer(f_fasta):
    RefFasta = dict()
    with open(f_fasta) as f:
        for line in f:
            if line.startswith("@"):
                id = line.strip().split("@")[1]
            else:
                seq = line.strip()
                RefFasta[id] = seq
    return RefFasta

def record_mapping(f_bam):
    Bam = pysam.AlignmentFile(f_bam, "rb")
    BamInfo = defaultdict(list)
    for read in Bam:
        read_name = read.qname
        chrom = read.reference_name
        cigar = read.cigartuples
        blocks = read.get_blocks()
        start = min([x[0] for x in blocks])
        end = max([x[1] for x in blocks])
        BamInfo[read_name].append([chrom, start, end, cigar, read])
    return BamInfo

def get_pos_pair(seq1, start1, seq2, start2, posStart, posEnd):
    seq1 = list(seq1)
    seq2 = list(seq2)

    posStart_to_start_dis = posStart-start1
    posStart_ref_num = seq1[:posStart_to_start_dis].count("-")
    posStart_query_num = seq2[:posStart_to_start_dis].count("-")
    rel_start = start2 + posStart_to_start_dis + posStart_ref_num - posStart_query_num

    posEnd_to_start_dis = posEnd-start1
    posEnd_ref_num = seq1[:posEnd_to_start_dis].count("-")
    posEnd_query_num = seq2[:posEnd_to_start_dis].count("-")
    rel_end = start2 + posStart_to_start_dis + posEnd_ref_num - posEnd_query_num

    return rel_start, rel_end

def sequence_pair(query_seq, ref_seq, cigartuple):
    query_pair_seq = list()
    ref_pair_seq = list()
    extrat_pair = 0
    extrat_ref = 0
    for cigar, length in cigartuple:
        if cigar == 0:
            p_start, p_end = extrat_pair, extrat_pair + length
            query_pair_seq.extend(query_seq[p_start:p_end])
            extrat_pair += length

            r_start, r_end = extrat_ref, extrat_ref+length
            ref_pair_seq.extend(ref_seq[r_start:r_end])
            extrat_ref += length

        elif cigar == 1:
            query_pair_seq.extend(["-" for i in range(length)])
            r_start, r_end = extrat_ref, extrat_ref+length
            ref_pair_seq.extend(ref_seq[r_start:r_end])
            extrat_ref += length
            
        elif cigar == 2:
            ref_pair_seq.extend(["-" for i in range(length)])
            p_start, p_end = extrat_pair, extrat_pair + length
            query_pair_seq.extend(query_seq[p_start:p_end])
            extrat_pair += length
        elif cigar == 4:
            query_pair_seq.extend(["-" for i in range(length)])
            r_start, r_end = extrat_ref, extrat_ref+length
            ref_pair_seq.extend(ref_seq[r_start:r_end])
            extrat_ref += length
        else:    
            print(read, cigar)

    return "".join(query_pair_seq), "".join(ref_pair_seq)


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-f", "--flank", type=str, dest="flank_fasta", metavar="cottonMD_site.fasta", required=True)
    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="Ghirsutum_HAU_genome.fasta", required=True)
    base_group.add_argument("-b", "--bam", type=str, dest="bam", metavar="GWAS.bam", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("-s", "--sequence", type=str, dest="sequence", metavar="GWAS_pair_seq.txt", required=True)
    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_fasta = args.flank_fasta
    f_genome = args.genome
    f_bam = args.bam
    f_out = args.output
    f_out_seq = args.sequence

    Fasta = pysam.FastaFile(f_genome)
    RefFasta = record_refer(f_fasta)
    BamInfo = record_mapping(f_bam)

    with open(f_out_seq, "w") as f_o1, open(f_out, "w") as f_o2:
        for read in BamInfo:
            ref_seq = RefFasta[read]
            for info in BamInfo[read]:
                chrom, start, end,  cigartuple, read_info = info
                query_seq = list(Fasta.fetch(chrom, start, end))
                query_pair_seq, ref_pair_seq = sequence_pair(query_seq, ref_seq, cigartuple)
                
                query_head = "%s:%s:%s:: %s" % (chrom, start, end, read_info.cigarstring)
                f_o1.write(read+"\n")
                f_o1.write(ref_pair_seq+"\n")
                f_o1.write(query_head+"\n")
                f_o1.write(query_pair_seq+"\n")

                ref_start, ref_end = read.split("::")[-1].split(":")
                ref_chrom, ref_flank_start, ref_flank_end = read.split("::")[1].split(":")
                rel_start, rel_end = get_pos_pair(ref_pair_seq, int(ref_flank_start), query_pair_seq, int(start), int(ref_start), int(ref_end))
                outline = [ref_chrom, ref_start, ref_end, chrom, rel_start, rel_end, read_info.cigarstring]
                outline = list(map(str, outline))
                f_o2.write("\t".join(outline)+"\n")

            
if __name__ == "__main__":
    main()




