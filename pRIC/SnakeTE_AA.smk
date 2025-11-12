
rule TE:
    input:
        "./results/TE/EP.pair.TE.statis.txt",


rule EP_pair_with_TE:
    input:
        EP = "./results/EP/A_WT_merge/A_WT_merge.EP.fisher.txt",
        TE_anno = "../Supplymental/Gabor_genome/TE/Ga.repeat.gff3",
    output:
        res = "./results/TE/EP.pair.TE.txt",
        statis = "./results/TE/EP.pair.TE.statis.txt",
    shell:
        """
        python scripts/TE/fetch_TE_in_EP_pair.py --EP {input.EP} -TE {input.TE_anno} -o {output.res}  -s {output.statis}
        """
       