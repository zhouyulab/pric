include: "MetaInfo_Tetraploid.smk"

rule mapping:
    input:
        expand("results/chimeric/{sample}/{sample}.chimeric.STAR.txt", sample=SAMPLEs),
        expand("results/chimeric/{sample}/{sample}.chimeric.BWA.txt", sample=SAMPLEs),
        expand("results/chimeric/{sample}/{sample}.chimeric.filter.txt", sample=SAMPLEs),

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
        python scripts/fetch_chimeric_from_STAR_aligned.py -b {input.bam} -c {input.chimeric} \
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
        python scripts/fetch_chimeric_from_BWA_aligned.py -b {input.bam} -R {input.fastq} \
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
        python scripts/merge_chimeric.py -i {output.merge} -o {output.txt}
        """
