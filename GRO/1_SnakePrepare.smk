import os
# include: "AD_MetaInfo.smk"
include: "AA_MetaInfo.smk"


rule Process:
    input:
        expand("./results/prepare/cutadapt/{sample}.R1.fastq.gz", sample=SAMPLEs),
        expand("./results/prepare/fastqc1/{sample}.R2_fastqc.html", sample=SAMPLEs),
        expand("./results/prepare/trim/{sample}/{sample}.R1_val_1.fq.gz", sample=SAMPLEs),
        expand("./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz", sample=SAMPLEs),


rule cutadapt:
    input:
        fq1 = os.path.join(data_indir, "{sample}.R1.fastq.gz"),
        fq2 = os.path.join(data_indir, "{sample}.R2.fastq.gz"),
    output:
        fq1 = "./results/prepare/cutadapt/{sample}.R1.fastq.gz",
        fq2 = "./results/prepare/cutadapt/{sample}.R2.fastq.gz",
    threads: 8
    shell:
        """
        cutadapt --trim-n --match-read-wildcards -u 16 -n 4 -a AGATCGGAAGAGCACACGTCTG -a AAAAAAAA -a GAACTCCAGTCAC -e 0.2 --nextseq-trim 20 -m 15 -j {threads} \
        -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2}					
        """
    
rule fastqc:
    input:
        fq1 = "./results/prepare/cutadapt/{sample}.R1.fastq.gz",
        fq2 = "./results/prepare/cutadapt/{sample}.R2.fastq.gz",
    output:
        html = "./results/prepare/fastqc1/{sample}.R2_fastqc.html"
    params:
        odir = "./results/prepare/fastqc1/"
    threads: 5
    log: "./results/prepare/fastqc1/{sample}.log",
    shell:"""
        fastqc -o {params.odir} -t {threads} {input.fq1} &>{log}
        fastqc -o {params.odir} -t {threads} {input.fq2} &>{log}
        """

rule trim:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
    output:
        fq1 = "./results/prepare/trim/{sample}/{sample}.R1_val_1.fq.gz",
        fq2 = "./results/prepare/trim/{sample}/{sample}.R2_val_2.fq.gz",
    params:
        prefix = "./results/prepare/trim/{sample}"
    threads: 8
    shell:
        """
        trim_galore --paired --quality 20 --length 15 -j {threads} --fastqc --small_rna -o {params.prefix} {input.fq1} {input.fq2}
        """

rule rm_rRNA:
    input:
        fastq_R1 = "./results/prepare/trim/{sample}/{sample}.R1_val_1.fq.gz",
        fastq_R2 = "./results/prepare/trim/{sample}/{sample}.R2_val_2.fq.gz",
    output:
        bam = "./results/prepare/rmrRNA/{sample}_mapping/Aligned.sortedByCoord.out.bam",
        unmap_R1 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
        unmap_R2 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
    threads: 30
    params:
        index = rRNA_INDX,
        mapping_dir = "./results/prepare/rmrRNA/{sample}_mapping"
    shell:"""
        if [[ -e {params.mapping_dir} ]]; then
            rm -r {params.mapping_dir}
        fi
        mkdir -p {params.mapping_dir}

        STAR --runThreadN {threads} \
        --genomeDir {params.index} \
        --outFileNamePrefix {params.mapping_dir}/ \
        --readFilesIn {input.fastq_R1} {input.fastq_R2} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 100 --outSAMattributes All \
        --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
        --chimSegmentMin 15 --chimJunctionOverhangMin 15 \
        --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outReadsUnmapped Fastx \
        --limitOutSJcollapsed 500000000 \
        --readFilesCommand gunzip -c --limitBAMsortRAM 23056060626
        
        gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
        gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
        rm {params.mapping_dir}/Unmapped.out.mate1 {params.mapping_dir}/Unmapped.out.mate2
    """
