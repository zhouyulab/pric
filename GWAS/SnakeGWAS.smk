
Flanks = [200, 500, 1000, 2000, 3000, 5000, 10000]
Samples = ["AD_fl_merge", "AD_WT_merge", "AD_pCp_merge"]
rule GWAS:
    input:
        "./results/GWAS.chimeric.txt"
        "results/GWAS.random.txt",
        expand("results/{sample}_EP_GWAS.txt", sample=Samples),
        "results/GWAS.edge.txt",
        "results/eQTL.bin.txt",
        "results/eQTL.random.bin.txt",
        "results/eQTL.plot.txt",


rule GWAS_chimeric:
    input:
        chimeric =  "../pRIC/results/chimeric/AD_pCp_merge/AD_pCp_merge.chimeric.cnt.txt",
        GWAS = "./results/SNP.merge.txt",
    output:
        txt = "./results/GWAS.chimeric.txt"
    shell:
        """
        python scripts/GWAS_zscore.py -c {input.chimeric} -g {input.GWAS} -o {output.txt} 
        """

rule GWAS_random:
    input:
        chimeric = "../pRIC/results/chimeric/AD_pCp_merge/AD_pCp_merge.chimeric.cnt.txt",
        GWAS =  "./results/SNP.merge.txt",
        chrom_size = "../Supplymental/Ghir_genome/Ghirsutum.chrom.sizes"
    output:
        txt = "results/GWAS.random.txt"
    shell:
        """
        python scripts/GWAS_random.py -c {input.chimeric} -chrom {input.chrom_size} -g {input.GWAS} -o {output.txt} 
        """

rule integrat_GWAS_with_EP:
    input:
        SNP = "./results/SNP.merge.txt",
        EP = "../pRIC/results/EP/{sample}/{sample}.EP.anno.txt"
    output:
        txt = "results/{sample}_EP_GWAS.txt"
    shell:
        """
        python scripts/interagted_SNP_with_EP.py {input.SNP} {input.EP} {output.txt}
        """

rule anno_GWAS:
    input:
        txt = "results/AD_pCp_merge_EP_GWAS.txt",
    output:
        edge = "results/GWAS.edge.txt",
        node = "results/GWAS.node.txt",
    shell:
        """
        python scripts/GWAS_info_to_PPI.py {input.txt} {output.edge} {output.node}
        """

rule fetch_eQTL_chimeric:
    input:
        bed = "../Supplymental/Ghir_genome/Ghir.gene.bed6",
        QTL = "./results/eQTL.txt",
        chimeric = "../pRIC/results/chimeric/AD_pCp_merge/AD_pCp_merge.chimeric.cnt.txt",
    output:
        txt = "results/eQTL.bin.txt"
    shell:
        """
        python scripts/eQTL_zscore.py -g {input.bed} -e {input.QTL} -c {input.chimeric} -o {output.txt} 
        """

rule fetch_random_eQTL_chimeric:
    input:
        GWAS = "./results/eQTL.random.txt",
        gene = "../Supplymental/Ghir_genome/Ghir.gene.bed6",
        chimeric = "../pRIC/src/results/chimeric/AD_pCp_merge/AD_pCp_merge.chimeric.cnt.txt",
    output:
        txt = "results/eQTL.random.bin.txt",
    shell:
        """
        python scripts/eQTL_zscore_random.py -e {input.GWAS} -g {input.gene} -c {input.chimeric} -o {output.txt}
        """

rule merge_eQTL_and_random:
    input:
        eQTL = "results/eQTL.bin.txt",
        random =  "results/eQTL.random.bin.txt",
    output:
        txt = "results/eQTL.plot.txt"
    shell:
        """
        python scripts/plot_eQTL_zcore.py {input.eQTL} {input.random} {output.txt}
        """

rule statis_GWAS:
    input:
        GWAS = "./results/SNP.merge.txt",
        gene = "../Supplymental/Ghir_genome/Ghirsutum_HAU_gene_model.bed",
    output:
        txt = "results/GWAS.pos.anno.txt"
    shell:
        """
        python scripts/GWAS_Pos.anno.py -gene {input.gene} \
        -GWAS {input.GWAS} -o {output.txt}
        """