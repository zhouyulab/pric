import os
include: "AD_MetaInfo.smk"
# include: "AA_MetaInfo.smk"

rule Mapping:
    input:
        expand("./results/peak/{sample}.meta.csv", sample=SAMPLEs),
        expand("./results/peak/{sample}_Forward_peaks.xls", sample=SAMPLEs),
        "./results/peak/AD.TSS.signal.pdf"
        "./results/peak/AD_GRO.gtf",
        "./results/peak/AD_GRO.filter.bed",
        "../results/peak/AA.TSS.signal.pdf",

rule meta:
    input:
        bw_F = "./results/bw/{sample}.Forward.bw",
        bw_R = "./results/bw/{sample}.Reverse.bw",
        bed = GENE,
        chrom_size = CHROM_SIZE,
    output:
        txt = "./results/peak/{sample}.meta.csv",
    params:
        after = 3000,
        before = 3000,
        around = 0,
        bin_num_up = 50,
        bin_num_body = 120,
        bin_num_down = 50
    shell:
        """
        python scripts/analysis/meta.py --input {input.bed} --bw-R {input.bw_R} --bw-F {input.bw_F} \
        --output {output.txt}  --after {params.after} --before {params.before} --chrom_size {input.chrom_size} \
        --around {params.around} --bin_num_up {params.bin_num_up} --bin_num_body {params.bin_num_body} --bin_num_down {params.bin_num_down}
        """

rule peak:
    input:
        bam_F = "./results/noNorm/{sample}.Forward.bam",
        bam_R = "./results/noNorm/{sample}.Reverse.bam",
    params:
        content = "./results/peak/",
        prefix_F = "{sample}_Forward",
        prefix_R = "{sample}_Reverse",
    output:
        peak_F = "./results/peak/{sample}_Forward_peaks.xls",
        peak_R = "./results/peak/{sample}_Reverse_peaks.xls",
    log: "./results/peak/{sample}.log",
    shell:
        """
        macs2 callpeak -t {input.bam_F} -f BAM --nomodel -p 0.05 -n {params.prefix_F} --keep-dup all --outdir {params.content} &> {log}
        macs2 callpeak -t {input.bam_R} -f BAM --nomodel -p 0.05 -n {params.prefix_R} --keep-dup all --outdir {params.content} &> {log}
        """

rule AD_TSS:
    input:
        bed = GENE,
        H3K27ac = "../AD_ChIP/results/bw/WT_H3K27ac.bw",
        K3K4me3 = "../AD_ChIP/results/bw/WT_H3K4me3.bw",
        DNase = "../DNase/results/bw/DNase_Rep1.bw",
        GRO = "./results/bw/GRO_AD_WT_Rep1.bw"
    output:
        mat = "./results/peak/AD.TSS.signal.txt.gz"
    threads: 8
    shell:
        """
        computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
        -R {input.bed} -S {input.H3K27ac} {input.K3K4me3} {input.DNase} {input.GRO} \
       --skipZeros -o {output.mat} --numberOfProcessors {threads}
        """

rule AD_TSS_plot:
    input:
        mat = "./results/peak/AD.TSS.signal.txt.gz"
    output:
        pdf = "./results/peak/AD.TSS.signal.pdf"
    shell:
        """
        plotProfile -m {input.mat} --perGroup -out {output.pdf} --colors red blue green orange --plotType lines --plotHeight 7 --plotWidth 9
        """

rule merge_GRO_AD:
    input:
        peak_rep1_F = "./results/peak/GRO_AD_WT_Rep1_Forward_peaks.narrowPeak",
        peak_rep1_R = "./results/peak/GRO_AD_WT_Rep1_Reverse_peaks.narrowPeak",
        peak_rep2_F = "./results/peak/GRO_AD_WT_Rep2_Forward_peaks.narrowPeak",
        peak_rep2_R = "./results/peak/GRO_AD_WT_Rep2_Reverse_peaks.narrowPeak",
    output:
        bed = "./results/peak/AD_GRO.bed",
        gtf = "./results/peak/AD_GRO.gtf",
    shell:
        """
        python scripts/analysis/merge_peak.py --prak-R-1 {input.peak_rep1_R} --peak-F-1 {input.peak_rep1_F} \
        --prak-R-2 {input.peak_rep2_R} --peak-F-2 {input.peak_rep2_F} -o {output.bed} -g {output.gtf}
        """

rule filter_GRO_peak_AD:
    input:
        gtf = "./results/peak/AD_GRO.gtf",
        bam = expand("./results/uniq/{sample}.uniq.bam", sample = ["GRO_AD_WT_Rep1", "GRO_AD_WT_Rep2"]),
    output:
        txt = "./results/peak/AD_GRO.featureCount.txt",
        summary = "./results/peak/AD_GRO.featureCount.txt.summary",
        bed = "./results/peak/AD_GRO.filter.bed",
    threads: 10
    shell:
        """
        featureCounts -s 1 -a {input.gtf} -o {output.txt}  -t exon -g gene_name --minOverlap 3 -O -p -C -T {threads} {input.bam}
        python scripts/analysis/filter_peak.py -i {output.txt} -s {output.summary} -o {output.bed}
        """

rule AA_TSS:
    input:
        bed = "../Supplymental/Gabor_genome/Gabor.bed"
        H3K27ac = "../AA_ChIP/results/bw/H3K27ac.bw",
        K3K4me3 = "../AA_ChIP/results/bw/H3K4me3.bw",
        DNase = "../DNase/results/bw/AA_DNase_Rep1.bw",
        GRO = "../results/bw/GRO_AA_WT_Rep2.bw"
    output:
        mat = "../results/peak/AA.TSS.signal.txt.gz"
    threads: 8
    shell:
        """
        computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
        -R {input.bed} -S {input.H3K27ac} {input.K3K4me3} {input.DNase} {input.GRO} \
       --skipZeros -o {output.mat} --numberOfProcessors {threads}
        """

rule AA_TSS_plot:
    input:
        mat = "../results/peak/AA.TSS.signal.txt.gz"
    output:
        pdf = "../results/peak/AA.TSS.signal.pdf"
    shell:
        """
        plotProfile -m {input.mat} --perGroup -out {output.pdf} --colors red blue green orange --plotType lines --plotHeight 7 --plotWidth 9
        """