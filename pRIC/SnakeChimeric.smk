include: "MetaInfo_Tetraploid.smk"
merge_Sample = ["AD_WT_merge", "AD_fl_merge", "AD_pCp_merge", "A_WT_merge"]
rule mapping:
    input:
        expand("results/chimeric/{sample}/{sample}.chimeric.STAR.txt", sample=SAMPLEs),
        expand("results/chimeric/{sample}/{sample}.chimeric.BWA.txt", sample=SAMPLEs),
        expand("results/chimeric/{sample}/{sample}.chimeric.filter.txt", sample=SAMPLEs),
        expand("../results/chimeric/{group}_merge/{group}_merge.chimeric.filter.txt", group=["A_WT", "AD_WT"]),
        expand("results/chimeric/{sample}/{sample}.chimeric.cnt.txt",sample=SAMPLEs+merge_Sample)

rule fetch_chimeric_from_STAR:
    input:
        bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.sort.uniq.bam",
        chimeric = "results/mapping/{sample}/STAR_mapping/{sample}.Chimeric.out.sam",
        bed = Gh_INTRON,
        fastq_R1 = "../results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
        fastq_R2 = "../results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
    params:
        rule = '1+-,1-+,2++,2--',
    output:
        txt = "results/chimeric/{sample}/{sample}.chimeric.STAR.txt",
    shell:
        """
        python scripts/chimeric/fetch_chimeric_from_STAR_aligned.py -b {input.bam} -c {input.chimeric} \
        -R1 {input.fastq_R1} -R2 {input.fastq_R2} --known-js {input.bed} --rule {params.rule} -o {output.txt}
        """

rule fetch_chimeric_from_BWA:
    input:
        bam = "results/mapping/{sample}/BWA_mapping/{sample}.sam",
        bed = Gh_INTRON,
        fastq = "results/mapping/{sample}/BWA_mapping/{sample}_collapse.fq.gz",
    output:
        txt = "results/chimeric/{sample}/{sample}.chimeric.BWA.txt",
    shell:
        """
        python scripts/chimeric/fetch_chimeric_from_BWA_aligned.py -b {input.bam} -R {input.fastq} \
        --known-js {input.bed} -o {output.txt}
        """

rule merge_chimeric:
    input:
        STAR = "results/chimeric/{sample}/{sample}.chimeric.STAR.txt",
        BWA = "results/chimeric/{sample}/{sample}.chimeric.BWA.txt",
    output:
        merge = "results/chimeric/{sample}/{sample}.chimeric.merge.txt",
        txt = "results/chimeric/{sample}/{sample}.chimeric.filter.txt",
    shell:
        """
        head -n 1 {input.STAR} > {output.merge}
        sed 1d {input.STAR} >> {output.merge}
        sed 1d {input.BWA} >> {output.merge} 
        python scripts/chimeric/merge_chimeric.py -i {output.merge} -o {output.txt}
        """

rule merge_replication:
    input:
        rep1 = "results/chimeric/{group}_Rep1/{group}_Rep1.chimeric.filter.txt",
        rep2 = "results/chimeric/{group}_Rep2/{group}_Rep2.chimeric.filter.txt",
    output:
        txt = "results/chimeric/{group}_merge/{group}_merge.chimeric.filter.txt",
    shell:
        """
        cat {input.rep1} > {output.txt}
        sed 1d {input.rep2} >> {output.txt}
        """

rule merge_wt_and_fl:
    input:
        WT = "results/chimeric/AD_WT_merge/AD_WT_merge.chimeric.filter.txt",
        fl = "results/chimeric/AD_fl_merge/AD_fl_merge.chimeric.filter.txt",
    output:
        merge = "results/chimeric/AD_pCp_merge/AD_pCp_merge.chimeric.filter.txt",
    shell:
        """
        cat {input.WT} > {output.merge}
        sed 1d {input.fl} >> {output.merge}
        """


rule cnt_chimeric:
    input:
        txt = "results/chimeric/{sample}/{sample}.chimeric.filter.txt",
    output:
        txt = "results/chimeric/{sample}/{sample}.chimeric.cnt.txt",
    shell:
        """
        python scripts/chimeric/cnt_chim_reads.py -i {input.txt} -o {output.txt}
        """

