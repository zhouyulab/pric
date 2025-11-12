import os
# include: "AA_MetaInfo.smk"
include: "AD_MetaInfo.smk"


outdir = "./results/"
rule mapping:
    input:
        expand(os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.bam"),sample=SAMPLEs),
        expand(os.path.join(outdir, "uniq/{sample}.uniq.bam"),sample=SAMPLEs),
        expand("./results/mapping/{sample}_mapping/{sample}.infer.txt",sample=SAMPLEs),
        expand(os.path.join(outdir, "bw/{sample}.bw"),sample=SAMPLEs),
        expand("./results/peak/{sample}_DNase_peaks.xls",sample=SAMPLEs),
        "./results/peak/AD_DNase_merge.peak",

rule STAR_mapping:
    input:
        fq1 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
        fq2 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
    params:
        index = STAR_INDX,
        prefix = os.path.join(outdir, "mapping/{sample}_mapping/"),
    threads: 20
    output:
        os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.bam"),
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
        bam = os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.bam"),
    output:
        bam = os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.sort.bam"),
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule uniq:
    input:
        bam = os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.sort.bam"),
    output:
        bam = os.path.join(outdir, "uniq/{sample}.uniq.bam"),
        txt = os.path.join(outdir, "uniq/{sample}_metrics.txt")
    log: os.path.join(outdir, "uniq/{sample}.log")
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true  --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.txt} &> {log}
        samtools index {output.bam}
        """

rule bam2wig:
    input:
        bam = rules.uniq.output.bam,
        chromSize = CHROM_SIZE
    output:
        bw = os.path.join(outdir, "bw/{sample}.bw"),
    threads:1
    params:
        prefix = os.path.join(outdir, "bw/{sample}"),
    shell:
        """
        bam2wig.py -i {input.bam} -s {input.chromSize} -u -t 100000000 -o {params.prefix}
        """

rule infer_experiment:
    input:
        bam = os.path.join(outdir, "mapping/{sample}_mapping/Aligned.sortedByCoord.out.bam"),
        Bed = BED,
    output:
        txt = "./results/mapping/{sample}_mapping/{sample}.infer.txt",
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.Bed} > {output.txt}
        """

rule macs2:
    input:
        bam = rules.uniq.output.bam,
    output:
        peak = "./results/peak/{sample}_DNase_peaks.xls",
    params:
        prefix = "{sample}_DNase",
        content = "./results/peak/",
    log:
        "./results/peak/{sample}.peak.log"
    shell:
        """
        macs2 callpeak -t {input.bam} -f BAM --nomodel -p 0.05 -n {params.prefix} --outdir {params.content} &> {log}
        """

rule merge_peak:
    input:
        rep1 = lambda wildcards: "./results/peak/%s_DNase_peaks.narrowPeak" % SAMPLEs[0],
        rep2 = lambda wildcards: "./results/peak/%s_DNase_peaks.narrowPeak" % SAMPLEs[1],
    output:
        peak = "./results/peak/AD_DNase_merge.peak",
    shell:
        """
        python scripts/peaks/merge_peak.py --input-rep1 {input.rep1} --input-rep2 {input.rep2} -o {output.peak}
        """