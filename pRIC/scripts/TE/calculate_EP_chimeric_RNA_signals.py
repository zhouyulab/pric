import pyBigWig
import sys, argparse
import numpy as np

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.txt", required=True)
    base_group.add_argument("-b", "--bw", type=str, dest="bw", metavar="sample.bw", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    return(parser.parse_args(args))

def fetch_signal(chrom, site, bw, flank):
    site = int(site)
    start = max(0, site-flank)
    end = site + flank
    bw_values = bw.stats(chrom, start, end, type="mean")
    bw_values = np.nan_to_num(bw_values)
    assert len(bw_values) == 1
    if bw_values[0] is None:
        return 0
    else:
        return bw_values[0]

def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_input = args.input
    f_out = args.output
    f_bw = args.bw

    bw = pyBigWig.open(f_bw)
    flank = 5
    with open(f_input) as fi, open(f_out, "w") as fo:
        for line in fi:
            data = line.strip().split("\t")
            if data[0] == "DonorChrom":
                header = data + ["PromoterValue", "EnhancerValue"]
                fo.write("\t".join(header)+"\n")
                continue
            DonorChrom,DonorCJS,DonorStrand,DonorBlocks,AcceptorChrom,AcceptorCJS,AcceptorStrand,AcceptorBlocks,ReadNum,DonorFeature, AcceptorFeature,Group = data
            if Group != "E-P":
                continue
            if DonorFeature.startswith("pro"):
                promoterChrom, promoterCJS, promoterStrand, promoterBlocks = DonorChrom,DonorCJS,DonorStrand,DonorBlocks
                enhancerChrom, enhancerCJS, enhancerStrand, enhancerBlocks = AcceptorChrom,AcceptorCJS,AcceptorStrand,AcceptorBlocks
            else:
                assert AcceptorFeature.startswith("pro")
                promoterChrom, promoterCJS, promoterStrand, promoterBlocks = AcceptorChrom,AcceptorCJS,AcceptorStrand,AcceptorBlocks
                enhancerChrom, enhancerCJS, enhancerStrand, enhancerBlocks = DonorChrom,DonorCJS,DonorStrand,DonorBlocks
            promoter_value = fetch_signal(promoterChrom, promoterCJS, bw, flank)
            enhancer_value = fetch_signal(enhancerChrom, enhancerCJS, bw, flank)

            outline = [promoterChrom, promoterCJS, promoterStrand, promoterBlocks, enhancerChrom, enhancerCJS, enhancerStrand, enhancerBlocks,ReadNum,DonorFeature, AcceptorFeature,Group ]
            outline = outline + [str(promoter_value), str(enhancer_value)]
            fo.write("\t".join(outline)+"\n")

if __name__ == "__main__":
    main()
