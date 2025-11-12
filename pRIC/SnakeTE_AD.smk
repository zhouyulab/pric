include: "MetaInfo_Tetraploid.smk"
rule TE:
    input:
        expand("./results/TE/{sample}/PRE_EP.pair.TE.statis.txt", sample=["AD_pCp_merge"]),
        expand("./results/TE/chimeric/EP.chimeric.{flank}.TE.txt", flank=[100, 200, 250, 500, 1000]),
        expand("./results/TE/chimeric/{sample}_{rep}_random_EP.chimeric.{flank}.TE.txt", sample=["WT", "fl"], rep=["Rep1", "Rep2"], flank=[250]),
        expand("./results/TE/chimeric/{sample}_{rep}_{flank}_random_EP.MFE.txt", sample=["WT", "fl"], rep=["Rep1", "Rep2"], flank=[250]),


rule EP_pair_with_TE_part:
    input:
        EP = "./results/EP/{sample}/{sample}.EP.fisher.txt",
        TE_anno = "../Supplymental/Ghir_genome/TE/genome/genome.fa.out.tsv",
    output:
        res = "./results/TE/{sample}/PRE_EP.pair.TE.txt",
        statis = "./results/TE/{sample}/PRE_EP.pair.TE.statis.txt",
    shell:
        """
        python scripts/TE/fetch_TE_in_EP_pair.py --EP {input.EP} -TE {input.TE_anno} -o {output.res}  -s {output.statis}
        """

rule EP_chimeric_with_TE:
    input:
        EP_chimeric = "./results/EP/AD_pCp_merge/AD_pCp_merge.EP.chimeric.txt",
        TE_anno = "../Supplymental/Ghir_genome/TE/genome/genome.fa.out.tsv",
    params:
        flank = '{flank}'
    output:
        res =  "./results/TE/chimeric/EP.chimeric.{flank}.TE.txt",
        statis =  "./results/TE/chimeric/EP.chimeric.{flank}.TE.statis.txt",
    shell:
        """
        python scripts/TE/fetch_TE_in_EP_chimeric.py -TE {input.TE_anno} -c {input.EP_chimeric} \
        --flank {params.flank} -o {output.res} -s {output.statis}
        """

rule EP_chimeric_with_TE_random:
    input:
        bw = "../RNA-seq/results/bw/{sample}_{rep}.bw",
        chrom = CHROM_SIZE,
        TE_anno = "../Supplymental/Ghir_genome/TE/genome/genome.fa.out.tsv",
    params:
        flank = '{flank}'
    output:
        random_EP = "./results/TE/chimeric/{sample}_{rep}_{flank}_random_EP.chimeric.txt",
        res =  "./results/TE/chimeric/{sample}_{rep}_random_EP.chimeric.{flank}.TE.txt",
        statis =  "./results/TE/chimeric/{sample}_{rep}_random_EP.chimeric.{flank}.TE.statis.txt",
    shell:
        """
        python scripts/TE/random_EP_chimeric.py -c {input.chrom} -b {input.bw} -o {output.random_EP} -r 10000
        python scripts/TE/fetch_TE_in_EP_chimeric.py -TE {input.TE_anno} -c  {output.random_EP} \
        --flank {params.flank} -o {output.res} -s {output.statis}
        """

rule random_MFE:
    input:
        txt = "./results/TE/chimeric/{sample}_{rep}_{flank}_random_EP.chimeric.txt",
        fasta = fasta,
    params:
        temp_dir = "./results/TE/temp/temp_{sample}_{rep}_{flank}",
        flank = "{flank}"
    output:
        txt = "./results/TE/chimeric/{sample}_{rep}_{flank}_random_EP.MFE.txt",
    shell:
        """
        mkdir {params.temp_dir}
        python scripts/TE/calculate_random_EP_MFE.py -g {input.fasta} -t {params.temp_dir} -c {input.txt} \
        --flank {params.flank} -o {output.txt}
        """



