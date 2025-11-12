import os
include: "AA_MetaInfo.smk"

outdir = "./results/"
rule mapping:
    input:
        expand(os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.bam"),sample=SAMPLEs, ty=Types, rep=Reps),
        expand(os.path.join(outdir, "uniq/{sample}_{ty}_{rep}.uniq.bam"),sample=SAMPLEs, ty=Types, rep=Reps),
        expand("./results/mapping/{sample}_{ty}_{rep}_mapping/{sample}_{ty}_{rep}.infer.txt",sample=SAMPLEs, ty=Types, rep=Reps),
       

rule STAR_mapping:
    input:
        fq1 = os.path.join(data_indir, "{sample}_{ty}_{rep}_R1.fq.gz"),
        fq2 = os.path.join(data_indir, "{sample}_{ty}_{rep}_R2.fq.gz"),
    params:
        index = STAR_INDX,
        prefix = os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/"),
    threads: 20
    output:
        os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.bam"),
    shell:"""
        STAR --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesCommand zcat \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --alignEndsType Local \
        --outFilterType BySJout \
        --readFilesIn {input.fq1} {input.fq2} \
        --outFileNamePrefix {params.prefix} \
        --outStd Log --outSAMmode Full \
        --outFilterMismatchNmax 2 \
        --outFilterMultimapNmax 1 \
        --outFilterMultimapScoreRange 1 \
        --outSAMattributes All --outReadsUnmapped Fastx
        """

rule sort_bam:
    input:
        bam = os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.bam"),
    output:
        bam = os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.sort.bam"),
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule uniq:
    input:
        bam = os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.sort.bam"),
    output:
        bam = os.path.join(outdir, "uniq/{sample}_{ty}_{rep}.uniq.bam"),
        txt = os.path.join(outdir, "uniq/{sample}_{ty}_{rep}_metrics.txt")
    log: os.path.join(outdir, "uniq/{sample}_{ty}_{rep}.log")
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true  --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.txt} &> {log}
        samtools index {output.bam}
        """

rule infer_experiment:
    input:
        bam = os.path.join(outdir, "mapping/{sample}_{ty}_{rep}_mapping/Aligned.sortedByCoord.out.bam"),
        bed = Gh_BED,
    output:
        txt = "./results/mapping/{sample}_{ty}_{rep}_mapping/{sample}_{ty}_{rep}.infer.txt"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}
        """