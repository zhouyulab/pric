
rule conserved_enhancer:
    input:
        "./results/conserved/promoter_index/promoter.not",
        "./results/conserved/enhancer_index/enhancer.not",
        "./results/conserved/promoter_blast.out",
        "./results/conserved/enhancer_blast.out",
        "./results/conserved/promoter.pair.txt",
        "./results/conserved/enhancer.pair.txt",
        "./results/conserved/EP.conserved.pair.txt",

rule blast_index:
    input:
        enhancer = "./results/feature/enhancer.group.fa",
        promoter = "./results/feature/promoter.fa"
    output:
        promoter_index = "./results/conserved/promoter_index/promoter.not",
        enhancer_index = "./results/conserved/enhancer_index/enhancer.not",
    params:
        promoter_prefix = "./results/conserved/promoter_index/promoter",
        enhancer_prefix = "./results/conserved/enhancer_index/enhancer",
    shell:
        """
        makeblastdb -in {input.promoter} -dbtype nucl -parse_seqids -out {params.promoter_prefix}
        makeblastdb -in {input.enhancer} -dbtype nucl -parse_seqids -out {params.enhancer_prefix}
        """

rule blast:
    input:
        promoter_fasta = "../AD_ChIP/results/feature/promoter.fa",
        enhancer_fasta = "../AD_ChIP/results/feature/enhancer.group.fa",
    params:
        promoter_prefix = "./results/conserved/promoter_index/promoter",
        enhancer_prefix = "./results/conserved/enhancer_index/enhancer",
    output:
        promoter_out = "./results/conserved/promoter_blast.out",
        enhancer_out = "./results/conserved/enhancer_blast.out",
    threads: 6
    shell:
        """
        blastn -query {input.promoter_fasta} -db {params.promoter_prefix} -evalue 1e-6 -outfmt 6 -num_threads 6 -out {output.promoter_out}
        blastn -query {input.enhancer_fasta} -db {params.enhancer_prefix} -evalue 1e-6 -outfmt 6 -num_threads 6 -out {output.enhancer_out}
        """

rule conserved_promoter_and_enhancer:
    input:
        promoter_out = "./results/conserved/promoter_blast.out",
        enhancer_out = "./results/conserved/enhancer_blast.out",
    output:
        promoter = "./results/conserved/promoter.pair.txt",
        enhancer = "./results/conserved/enhancer.pair.txt",
    shell:
        """
        python scripts/conserved/fetch_conserved_promoter_enhancer.py -p {input.promoter_out} -e {input.enhancer_out} \
        --out-promoter {output.promoter} --out-enhancer {output.enhancer}
        """

rule EP_pair:
    input:
        promoter_pair = "./results/conserved/promoter.pair.txt",
        enhancer_pair = "./results/conserved/enhancer.pair.txt",
        A_EP = "../pRIC/AA_results/EP/A_WT_merge/A_WT_merge.EP.fisher.txt",
        AD_EP = "../pRIC/AD_results/EP/AD_pCp_merge/AD_pCp_merge.EP.fisher.txt"
    output:
        txt = "./results/conserved/EP.conserved.pair.txt",
    shell:
        """
        python scripts/conserved/fetch_conserved_EP.py -p {input.promoter_pair} -e {input.enhancer_pair} \
        --diploid {input.A_EP} --tetraploid {input.AD_EP} -o {output.txt}
        """