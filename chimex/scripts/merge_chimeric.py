import argparse
import sys
import logging
from collections import defaultdict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_input = args.input
    f_output = args.output

    chimeric_dict = defaultdict(int)
    with open(f_input) as fi:
        for line in fi:
            data = line.strip().split("\t")
            if data[0] == "DonorChrom":
                continue
            donor_chrom, donor_site, donor_strand, acceptor_chrom, acceptor_site, acceptor_strand = data[0], int(data[1]), data[2], data[5], int(data[6]), data[7]
            donor_block, acceptor_block = data[3], data[8]
            read_name, read_source = data[-2], data[-1]
            key = (donor_chrom, donor_site, donor_strand,donor_block, acceptor_chrom, acceptor_site, acceptor_strand, acceptor_block, read_name)
            chimeric_dict[key] += 1

    header = ["DonorChrom", "DonorCJS", "DonorStrand", "DonorBlock", "AcceptorChrom", "AcceptorCJS", "AcceptorStrand", "AcceptorBlock", "ReadName"]
    with open(f_output, "w") as fo:
        fo.write("\t".join(header)+"\n")
        for key in chimeric_dict:
            read_nums = chimeric_dict[key]
            read_name = key[-1]
            if read_nums >1:
                logging.info(read_name+" repeatedly extracted the junction site")
            outline = list(map(str, key))
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
