import os
include: "AD_MetaInfo.smk"
# include: "AA_MetaInfo.smk"

rule Mapping:
    input:
        expand("./results/mapping/{sample}_STAR/{sample}.Aligned.out.bam", sample=SAMPLEs),
        expand("./results/mapping/{sample}_STAR/{sample}.Aligned.out.uniq.bam", sample=SAMPLEs),
        expand("./results/uniq/{sample}.uniq.bam", sample=SAMPLEs),
        expand("./results/mapping/{sample}_STAR/{sample}.infer.txt", sample=SAMPLEs),
        expand("./results/bw/{sample}.bw", sample=SAMPLEs),
        expand("./results/bw/{sample}.Forward.bw",sample=SAMPLEs),
        expand("./results/noNorm/{sample}.Forward.bam",sample=SAMPLEs),
        "./results/bw/PRO_AD.pdf", 

rule STAR:
    input:
        R1 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
        R2 = "./results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
    output:
        bam = "./results/mapping/{sample}_STAR/{sample}.Aligned.out.bam",
    params:
        index = STAR_INDX,
        gtf = GTF,
        outdir = "./results/mapping/{sample}_STAR/{sample}.",
    threads: 16
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --sjdbGTFfile {params.gtf} \
        --readFilesIn {input.R1} {input.R2} \
        --outFileNamePrefix {params.outdir} \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped Fastx \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMultimapNmax 100 \
        --outSAMattributes All \
        --readFilesCommand zcat \
        --limitOutSJcollapsed 50000000
        """


rule infer_experiment:
    input:
        bam = "./results/mapping/{sample}_STAR/{sample}.Aligned.out.bam",
        bed = BED,
    output:
        txt = "./results/mapping/{sample}_STAR/{sample}.infer.txt"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}
        """

rule sort:
    input:
        bam = "./results/mapping/{sample}_STAR/{sample}.Aligned.out.bam",
    output:
        bam = "./results/mapping/{sample}_STAR/{sample}.Aligned.out.uniq.bam",
    threads:
        10
    shell:
        """
        samtools view -@ {threads} -hb -q 30 -F 256 {input.bam} | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule rmdup:
    input:
        bam = "./results/mapping/{sample}_STAR/{sample}.Aligned.out.uniq.bam",
    output:
        bam = "./results/uniq/{sample}.uniq.bam",
        txt = "./results/uniq/{sample}_metrics.txt",
    log: "./results/uniq/{sample}.log"
    shell:
        """
        picard MarkDuplicates -Xmx100g --REMOVE_DUPLICATES true  --INPUT {input.bam} --OUTPUT {output.bam} -MAX_FILE_HANDLES 100 --METRICS_FILE {output.txt} &> {log}
        samtools index {output.bam}
        """

rule bam2wig:
    input:
        bam = rules.rmdup.output.bam,
        chromSize = CHROM_SIZE,
    output:
        bw = "./results/bw/{sample}.bw",
    threads:1
    params:
        prefix = "./results/bw/{sample}",
    shell:
        """
        bam2wig.py -i {input.bam} -s {input.chromSize} -u -t 100000000 -o {params.prefix}
        """

rule bam2wig_strand:
    input:
        bam = rules.rmdup.output.bam,
        chromSize = CHROM_SIZE,
    output:
        bw_R = "./results/bw/{sample}.Forward.bw",
        bw_F = "./results/bw/{sample}.Reverse.bw",
    threads:1
    params:
        prefix = "./results/bw/{sample}",
    shell:
        """
        bam2wig.py -i {input.bam} -s {input.chromSize} -d 1++,1--,2+-,2-+ -u -t 100000000 -o {params.prefix}
        """

rule split_bam:
    input:
        bam = rules.rmdup.output.bam,
    output:
        bw_F = "./results/noNorm/{sample}.Forward.bw",
        bw_R = "./results/noNorm/{sample}.Reverse.bw",
        bam_F = "./results/noNorm/{sample}.Forward.bam",
        bam_R = "./results/noNorm/{sample}.Reverse.bam",
    threads:10
    shell:
        """
        samtools view -@ {threads} -b -f 128 -F 16 {input.bam} > {wildcards.sample}.fwd1.bam
        samtools view -@ {threads} -b -f 80  {input.bam}  > {wildcards.sample}.fwd2.bam
        samtools merge -@ {threads} -f {output.bam_F} {wildcards.sample}.fwd1.bam {wildcards.sample}.fwd2.bam
        samtools index -@ {threads} {output.bam_F}
        bamCoverage --bam {output.bam_F} -o {output.bw_F}
        rm {wildcards.sample}.fwd1.bam {wildcards.sample}.fwd2.bam 

        samtools view -@ {threads} -b -f 144 {input.bam} > {wildcards.sample}.rev1.bam
        samtools view -@ {threads} -b -f 64 -F 16 {input.bam} > {wildcards.sample}.rev2.bam
        samtools merge -@ {threads} -f {output.bam_R} {wildcards.sample}.rev1.bam {wildcards.sample}.rev2.bam
        samtools index -@ {threads} {output.bam_R}
        bamCoverage --bam {output.bam_R} -o {output.bw_R}
        rm {wildcards.sample}.rev1.bam {wildcards.sample}.rev2.bam
        """

rule deeptools_rep:
    input:
        bw = expand("./results/bw/{sample}.bw", sample = SAMPLEs)
    output:
        npz = "./results/bw/GRO_AD.npz", 
    threads:
        10
    shell:
        """
        multiBigwigSummary bins -b {input} -o {output.npz} -p {threads}
        """

rule plot_replation:
    input:
        npz = "./results/bw/GRO_AD.npz", 
    output:
        pdf = "./results/bw/GRO_AD.pdf", 
        tab = "./results/bw/GRO_AD.tab", 
    shell:
        """
        plotCorrelation -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of ChIP" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        --removeOutliers --log1p \
        -o {output.pdf}  --outFileCorMatrix {output.tab}
        """

