import os
include: "AD_MetaInfo.smk"

rule mapping:
    input:
        expand("./results/bw/{sample}_{ty}.bw",sample=SAMPLEs, ty=Types),
        expand("./results/peak/{sample}_H3K4me3_peaks.narrowPeak",sample=SAMPLEs),
        expand("./results/peak/{sample}_H3K27ac_peaks.narrowPeak",sample=SAMPLEs),

rule merge_bam:
    input:
        bam1 = "./results/uniq/{sample}_{ty}_Rep1.uniq.bam",
        bam2 = "./results/uniq/{sample}_{ty}_Rep2.uniq.bam",
    output:
        bam = "./results/uniq/{sample}_{ty}_merge.bam",
    shell:
        """
        samtools merge {output.bam} {input.bam1} {input.bam2}
        samtools index {output.bam}
        """

rule bam2wig:
    input:
        bam = rules.merge_bam.output.bam,
        chromSize = CHROM_SIZE,
    output:
        bw = "./results/bw/{sample}_{ty}.bw",
    threads:10
    params:
        prefix = "./results/bw/{sample}_{ty}",
    shell:
        """
        bam2wig.py -i {input.bam} -s {input.chromSize} -u -t 100000000 -o {params.prefix}
        """

rule call_H3K4me3_peak:
    input:
        control = "./results/uniq/{sample}_Input_merge.bam",
        treatment = "./results/uniq/{sample}_H3K4me3_merge.bam",
    log: "./results/peak/{sample}_H3K4me3.log",
    output:
        xls = "./results/peak/{sample}_H3K4me3_peaks.xls",
        peak = "./results/peak/{sample}_H3K4me3_peaks.narrowPeak",
    params:
        perfix = "{sample}_H3K4me3",
        content =  "./results/peak/",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -f BAM --nomodel -q 0.01 -n {params.perfix} --outdir {params.content} &> {log}
        """

rule call_H3K27ac_peak:
    input:
        control = "./results/uniq/{sample}_Input_merge.bam",
        treatment = "./results/uniq/{sample}_H3K27ac_merge.bam",
    log: "./results/peak/{sample}_H3K27ac.log",
    output:
        xls = "./results/peak/{sample}_H3K27ac_peaks.xls",
        peak = "./results/peak/{sample}_H3K27ac_peaks.narrowPeak",
    params:
        perfix = "{sample}_H3K27ac",
        content =  "./results/peak/",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -f BAM --nomodel -q 0.01 -n {params.perfix} --outdir {params.content} &> {log}
        """

