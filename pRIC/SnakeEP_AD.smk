include: "MetaInfo_Tetraploid.smk"

rule EP:
    input:
        expand("./results/EP/{sample}/{sample}.EP.network", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/{sample}.EP.fisher.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/{sample}.EP.anno.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/{sample}.EP.nearest.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/{sample}.edge.txt",sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/enhancer.fa", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/streme_enhancer/streme.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/streme_promoter/streme.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),
        expand("./results/EP/{sample}/EP_flank_MFE.txt", sample=["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]),

rule fetch_EP_chimeric:
    input:
        chimeric = "./results/chmeric/{sample}/{sample}.chimeric.cnt.txt",
        enhancer = "../ChIP-seq/results/feature/enhancer.group.bed",
        promoter = "../ChIP-seq/results/feature/promoter.bed",
    output:
        txt = "./results/EP/{sample}/{sample}.EP.chimeric.txt",
        network = "./results/EP/{sample}/{sample}.EP.network",
    shell:
        """
        python scripts/EP/fetch_EP_chimeric.py -p {input.promoter} -e {input.enhancer} -c {input.chimeric} \
        -o {output.txt} -n {output.network}
        """

rule EP_fisher_test:
    input:
        network = rules.fetch_EP_chimeric.output.network
    output:
        txt = "./results/EP/{sample}/{sample}.all.fisher.txt",
        EP_txt = "./results/EP/{sample}/{sample}.EP.fisher.txt",
    params:
        cutoff = 2
    shell:
        """
        python scripts/EP/fisher.test.py -i {input.network} -o {output.txt} -c {params.cutoff} --out-EP {output.EP_txt}
        """

rule EP_anno:
    input:
        EP = "./results/EP/{sample}/{sample}.EP.fisher.txt",
        anno = "../Supplymental/Ghir_genome/gene.anno.fiber.tsv",
    output:
        txt = "./results/EP/{sample}/{sample}.EP.anno.txt",
    shell:
        """
        python scripts/EP/anno_EP_gene.py {input.EP} {input.anno} {output.txt}
        """

rule EP_PPI:
    input:
        chrom_size = CHROM_SIZE,
        EP = "./results/EP/{sample}/{sample}.EP.fisher.txt",
    output:
        edge = "./results/EP/{sample}/{sample}.edge.txt",
        node = "./results/EP/{sample}/{sample}.node.txt",
    shell:
        """
        python scripts/EP/EP_PPI.py -i {input.EP} -e {output.edge} -n {output.node} --chrom-size {input.chrom_size}
        """


rule EP_interaction_is_nearest:
    input:
        txt = "./results/EP/{sample}/{sample}.EP.fisher.txt",
        nearest = "../ChIP-seq/results/feature/enhancer_nearest_promoter.txt"
    output:
        txt = "./results/EP/{sample}/{sample}.EP.nearest.txt",
    shell:
        """
        python scripts/EP/judge_EP_nearest.py -i {input.txt} -a {input.nearest} -o {output.txt}
        """

rule fetch_EP_fasta:
    input:
        txt = "./results/EP/{sample}/{sample}.EP.chimeric.txt",
        fasta = fasta,
    output:
        promoter = "./results/EP/{sample}/promoter.fa",
        enhancer = "./results/EP/{sample}/enhancer.fa",
    params:
        prefix = "./results/EP/{sample}/",
        flank = 250
    shell:
        """
        python scripts/EP/fetch_EP_fasta.py -c {input.txt} -g {input.fasta} -pre {params.prefix} --flank {params.flank}
        """

rule promoter_Motif:
    input:
        fasta = rules.fetch_EP_fasta.output.promoter
    output:
        txt =  "./results/EP/{sample}/streme_promoter/streme.txt"
    params:    
        prefix = "./results/EP/{sample}/streme_promoter/"
    shell:
        """
        if [[ -e {params.prefix} ]]; then
            rm -rf {params.prefix}
        fi
        mkdir {params.prefix}
        streme -rna -nmotifs 10 --evalue 0.05 -minw 4 -maxw 30 --oc {params.prefix} --p {input.fasta}
        """

rule enhancer_Motif:
    input:
        fasta = rules.fetch_EP_fasta.output.enhancer
    output:
        txt =  "./results/EP/{sample}/streme_enhancer/streme.txt"
    params:    
        prefix = "./results/EP/{sample}/streme_enhancer/"
    shell:
        """
        if [[ -e {params.prefix} ]]; then
            rm -rf {params.prefix}
        fi
        mkdir {params.prefix}
        streme -rna -nmotifs 10 --evalue 0.05 -minw 4 -maxw 30 --oc {params.prefix} --p {input.fasta}
        """

rule calculate_MFE:
    input:
        fasta = fasta,
        chimeric = "./results/EP/{sample}/{sample}.EP.chimeric.txt",
    output:
        txt = "./results/EP/{sample}/EP_flank_MFE.txt"
    params:
        outdir = "./results/EP/{sample}/flank_MFE/",
        flank = 250
    shell:
        """
        if [[ -e {params.outdir} ]]; then
        rm -rf {params.outdir}
        fi
        mkdir {params.outdir}
        python scripts/EP/calculate_EP_MFE.py -g {input.fasta} -t {params.outdir} -c {input.chimeric} --flank {params.flank} -o {output.txt}
        """
