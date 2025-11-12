import pysam
import re
import os
import argparse
import sys
import gzip
from collections import defaultdict
import numpy as np

def load_sam(f_bam):
    ChimericSam = pysam.AlignmentFile(f_bam, "r")
    array = []
    for read in ChimericSam:
        if not array: 
            array.append(read)
        elif read.qname == array[0].qname:
            array.append(read)
        else:
            yield array
            array = [read]

def load_fastq(f_in):
    with gzip.open(f_in, "rb") as f_i:
        array = []
        for i, line_data in enumerate(f_i):
            line = line_data.rstrip().decode()
            array.append(line)
            if i % 4 == 3:
                yield array
                array = []

def load_chimeric_fastq(f_fastq):
    loader1 = load_fastq(f_fastq)
    read_fastq_dict = dict()
    while True:
        try:
            name, seq, strand, qua = next(loader1)
            name = name.strip().split(" ")[0].split("@")[1]
            read_fastq_dict[name] = seq
        except StopIteration:
            break
    return read_fastq_dict

def load_known_js(f_known_js):
    known_JS_dict = defaultdict(set)
    if f_known_js is not None:
        print("Loading known JS ...")
        with open(f_known_js) as fi:
            for line in fi:
                data = line.rstrip("\n").split("\t")
                chrom = data[0]
                start = int(data[1])
                end = int(data[2])
                known_JS_dict[chrom].add((start, end))
            all_js_num = sum([len(x) for x in known_JS_dict.values()])
        print("{0} known JS were loaded.".format(all_js_num))
    return known_JS_dict

def get_reverse_seq(sequence):
    base_dict = {'n':'N', "a":"T", "t":"A", "c":"G", "g":"C"}
    seq = list(sequence.lower())[::-1]
    rev_seq = [base_dict[base] for base in seq]
    return "".join(rev_seq)

def get_sequence_relative_pos(sequence, subseq):
    for x in re.finditer(subseq, sequence):
        start, end = x.span()
        return start+1, end

def get_read_strand(read, rule):
    if rule == "++,--":
        read_strand = "-" if read.is_reverse else "+"
    elif rule == "+-,-+":
        read_strand = "+" if read.is_reverse else "-"
    return read_strand
    

def fetch_chimeric_info(read_li, read_fastq):
    part1 = read_li[0]
    part2 = read_li[1]
    part1_strand = get_read_strand(part1, rule = "+-,-+")
    part2_strand = get_read_strand(part2, rule = "+-,-+")
    part1_block = ["%s-%s" % (b[0], b[1]) for b in part1.get_blocks()]
    part2_block = ["%s-%s" % (b[0], b[1]) for b in part2.get_blocks()]
    part1_sequence = get_reverse_seq(part1.query_alignment_sequence) if part1.is_reverse else part1.query_alignment_sequence
    part2_sequence = get_reverse_seq(part2.query_alignment_sequence) if part2.is_reverse else part2.query_alignment_sequence

    part1_start, part1_end = get_sequence_relative_pos(read_fastq, part1_sequence)
    part2_start, part2_end = get_sequence_relative_pos(read_fastq, part2_sequence)

    if part1_end < part2_start:
        donor_chrom, donor_strand, donor_block, donor_sequence = part1.reference_name, part1_strand, ";".join(part1_block), part1_sequence
        donor_site = part1.reference_end if part1_strand == "+" else part1.reference_start
        acceptor_chrom, acceptor_strand, acceptor_block, acceptor_sequence = part2.reference_name, part2_strand, ";".join(part2_block), part2_sequence
        acceptor_site = part2.reference_start if part2_strand == "+" else part2.reference_end
        mid_sequence = read_fastq[part1_end:(part2_start-1)]
    else:
        donor_chrom, donor_strand, donor_block, donor_sequence = part2.reference_name, part2_strand, ";".join(part2_block), part2_sequence
        donor_site = part2.reference_end if part2_strand == "+" else part2.reference_start
        acceptor_chrom, acceptor_strand, acceptor_block, acceptor_sequence = part1.reference_name, part1_strand, ";".join(part1_block), part1_sequence
        acceptor_site = part1.reference_start if part1_strand == "+" else part1.reference_end
        mid_sequence = read_fastq[part2_end:(part1_start-1)]
    return [donor_chrom, donor_site, donor_strand, donor_block, donor_sequence, acceptor_chrom, acceptor_site, acceptor_strand, acceptor_block, acceptor_sequence, mid_sequence, part1.qname, "BWA"]


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    
    base_group.add_argument("-b", "--bam", type=str, dest="Bam", metavar="Aligned.sort.uniq.bam", required=True)
    base_group.add_argument("-R", "--read", type=str, dest="read", metavar="merge.fastq.gz", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.fq.gz", required=True)
    base_group.add_argument("--known-js", type=str, dest="bed", metavar="ref.intron.bed", required=True)
    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_bam = args.Bam
    f_read = args.read
    f_output = args.output
    f_known_js = args.bed

    known_JS_dict = load_known_js(f_known_js)
    header = ["DonorChrom", "DonorCJS", "DonorStrand", "DonorBlock", "DonorSequence", "AcceptorChrom", "AcceptorCJS", "AcceptorStrand", "AcceptorBlock", "AcceptorSequence", "MidSequence", "ReadName", "Source"]
    fo = open(f_output, "w")
    fo.write("\t".join(header)+"\n")

    ChimericLoder = load_sam(f_bam)
    ChimericFastq = load_chimeric_fastq(f_read)

    while True:
        try:
            read_li = next(ChimericLoder)
        except StopIteration:
            break
        if len(read_li) == 1:
            continue
        read_name = read_li[0].qname
        read_fastq = ChimericFastq[read_name]
        read_sig = np.zeros(len(read_fastq))
        for read in read_li:
            # print(read)
            read_match_seq = read.query_alignment_sequence
            if read.is_reverse:
                read_match_seq = get_reverse_seq(read_match_seq)
            # print(read_match_seq)
            for x in re.finditer(read_match_seq, read_fastq):
                start, end = x.span()
                read_sig[start:end] += 1
        if np.array(read_sig>1).any():
            continue
        else:
            if len(read_li) == 2:
                outline = fetch_chimeric_info(read_li, read_fastq)
                donor_chrom, donor_site, acceptor_chrom, acceptor_site = outline[0], outline[1], outline[5], outline[6]
                if  donor_chrom != acceptor_chrom:
                    outline = list(map(str, outline))
                    fo.write("\t".join(outline)+"\n")
                else:
                    if (donor_site, acceptor_site) not in known_JS_dict[donor_chrom]:
                        outline = list(map(str, outline))
                        fo.write("\t".join(outline)+"\n")

    fo.close()

if __name__ == "__main__":
    main()
