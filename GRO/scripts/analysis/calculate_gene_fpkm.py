import argparse,sys

def get_total_reads(f_summary):
    ToatlReads = 0
    with open(f_summary) as f:
        for line in f:
            if line.startswith("Status"):
                continue
            lines = line.strip().split("\t")
            ToatlReads = ToatlReads + int(lines[-1])
    return ToatlReads

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="gene.txt", required=True)
    base_group.add_argument("-s", "--summary", type=str, dest="summary", metavar="gene.txt.summary", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)

    return(parser.parse_args(args))

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_in = args.input
    f_summary = args.summary
    f_out = args.output

    header = ['Geneid', 'Chrom', 'Start', 'End', 'Strand', 'Length', 'Count', 'FPKM']
    handle = open(f_out, "w")
    handle.write("\t".join(header)+"\n")

    ToatlReads = get_total_reads(f_summary)
    with open(f_in) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith("Geneid"):
                continue
            lines = line.strip().split("\t")
            Geneid = lines[0]
            Chrom = lines[1].split(";")[0]
            Start = min(lines[2].split(";"))
            End = max(lines[3].split(";"))
            Strand = lines[4].split(";")[0]
            Length = int(lines[5])
            Count = int(lines[6])
            FPKM = Count * 1e9 / Length / ToatlReads
            outline = [Geneid, Chrom, Start, End, Strand, Length, Count, FPKM]
            outline = list(map(str, outline))
            handle.write("\t".join(outline)+"\n")
    handle.close()

if __name__ == "__main__":
    main()
