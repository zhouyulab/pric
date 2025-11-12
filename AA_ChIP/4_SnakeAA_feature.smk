import os
include: "AA_MetaInfo.smk"

rule mapping:
    input:
        "results/feature/enhancer.bed",
        "results/feature/promoter.bed",
        "results/feature/enhancer.count.bed",
        "results/feature/enhancer.ROSE.pdf",
        "results/feature/enhancer_nearest_promoter.txt",
        "results/feature/Enhancer.region.pdf",
        "results/feature/Promoter.region.pdf",
        "results/feature/enhancer.group.fa",


rule merge_H3K27ac_with_ATAC_DNase_GRO:
    input:
        H3K27ac = "results/peak/AA_H3K27ac_peaks.narrowPeak",
        DNase = "../DNase/results/peak/AA_DNase_merge.peak",
        H3K4me4 = "results/peak/H3K4me3_Promoter.peak",
        GRO = "../AA_GRO/results/feature/AA_GRO.filter.bed",
    params:
        prefix = "results/feature/"
    output:
        enhancer = "results/feature/enhancer.bed",
        promoter = "results/feature/promoter.bed",
    shell:
        """
        python scripts/feature/merge_peak_with_DNase_GRO.py --DNase {input.DNase} --H3K4me4 {input.H3K4me4} \
        --H3K27ac {input.H3K27ac} --GRO {input.GRO} -e {output.enhancer} -p {output.promoter} -prefix {params.prefix}
        """

rule enhancer_sigle:
    input:
        peak = "results/feature/enhancer.bed",
        bam_rep1 = "results/uniq/H3K27ac_merge.bam",
        DNase_rep1 = "../DNase/results/uniq/AA_DNase_Rep1.uniq.bam",
        DNase_rep2 = "../DNase/results/uniq/AA_DNase_Rep2.uniq.bam",
    output:
        bed = "results/feature/enhancer.count.bed"
    shell:
        """
        bedtools multicov -bed {input.peak} -bams {input.bam_rep1} {input.DNase_rep1} {input.DNase_rep2} > {output.bed}
        """

rule class_enhancer:
    input:
        bed = "results/feature/enhancer.count.bed"
    output:
        bed = "results/feature/enhancer.group.bed",
        pdf = "results/feature/enhancer.ROSE.pdf",
    shell:
        """
        Rscript scripts/feature/class_enhancer_to_TE_and_SE.R {input.bed} {output.bed} {output.pdf}
        """

rule find_enhancer_nearst_promoter:
    input:
        promoter = "results/feature/promoter.bed",
        enhancer = "results/feature/enhancer.group.bed",
    output:
        txt = "results/feature/enhancer_nearest_promoter.txt",
    shell:
        """
        python scripts/feature/find_enhancer_nearst_promoter.py -p {input.promoter} -e {input.enhancer} -o {output.txt}
        """

rule promoter_feature:
    input:
        promoter = "results/feature/promoter.bed",
        H3K27ac = "results/bw/H3K27ac.bw",
        H3K4me4 = "results/bw/H3K4me3.bw",
        DNase = "../DNase/results/bw/AA_DNase_Rep1.bw",
        GRO_forward = "../GRO/results/bw/GRO_AA_WT_Rep2.Forward.bw",
        GRO_reverse = "../GRO/results/bw/GRO_AA_WT_Rep2.Reverse.bw",
    output:
        txt = "results/feature/Promoter.region.gz",
        matrix = "results/feature/Promoter.outFileNameMatrix",
    threads: 8
    params:
        flank = 1000
    shell:
        """
        computeMatrix scale-regions --regionsFileName {input.promoter} --scoreFileName {input.H3K27ac} {input.H3K4me4} {input.DNase} {input.GRO_forward} {input.GRO_reverse} \
        --outFileNameMatrix {output.matrix} --beforeRegionStartLength {params.flank} --afterRegionStartLength  {params.flank} --regionBodyLength 1000 --missingDataAsZero \
        --startLabel TSS --endLabel TES  --skipZeros --numberOfProcessors {threads} --outFileName {output.txt}
        """

rule plot_promoter:
    input:
        txt = "results/feature/Promoter.region.gz",
    output:
        png = "results/feature/Promoter.region.pdf",
        txt = "results/feature/Promoter.region.Matrix.gz",
    shell:
        """
        plotHeatmap -m {input.txt} -out {output.png} --colorMap Blues Greens Oranges RdBu RdBu Reds \
        --zMin 0 0 0 -3 -3 --zMax 0.5 3 0.2 3 3 --heatmapHeight 15 --heatmapWidth 6 --outFileNameMatrix {output.txt}
        """

rule enhancer_feature:
    input:
        enhancer = "results/feature/enhancer.group.bed",
        H3K27ac = "results/bw/H3K27ac.bw",
        H3K4me4 = "results/bw/H3K4me3.bw",
        DNase = "../DNase/results/bw/AA_DNase_Rep1.bw",
        GRO_forward = "../GRO/results/bw/GRO_AA_WT_Rep2.Forward.bw",
        GRO_reverse = "../GRO/results/bw/GRO_AA_WT_Rep2.Reverse.bw",
    output:
        txt = "results/feature/Enhancer.region.gz",
        matrix = "results/feature/Enhancer.outFileNameMatrix",
    threads: 8
    params:
        flank = 1000
    shell:
        """
        computeMatrix scale-regions --regionsFileName {input.enhancer} --scoreFileName {input.H3K27ac} {input.H3K4me4} {input.DNase} {input.GRO_forward} {input.GRO_reverse} \
        --outFileNameMatrix {output.matrix} --beforeRegionStartLength {params.flank} --afterRegionStartLength  {params.flank} --regionBodyLength 1000 --missingDataAsZero \
        --startLabel TSS --endLabel TES  --skipZeros --numberOfProcessors {threads} --outFileName {output.txt}
        """

rule plot_enhancer:
    input:
        txt = "results/feature/Enhancer.region.gz",
    output:
        png = "results/feature/Enhancer.region.pdf",
        txt = "results/feature/Enhancer.region.Matrix.gz",
    shell:
        """
        plotHeatmap -m {input.txt} -out {output.png} --colorMap Blues Greens Oranges RdBu RdBu Reds \
        --zMin 0 0 0 -0.4 -0.4 --zMax 0.5 0.5 0.2 0.3 0.3 --heatmapHeight 15 --heatmapWidth 6 --outFileNameMatrix {output.txt}
        """

rule feature_fasta:
    input:
        promoter = "results/feature/promoter.bed",
        enhancer = "results/feature/enhancer.group.bed",
        genome = fasta,
    output:
        promoter = "results/feature/promoter.fa",
        enhancer = "results/feature/enhancer.group.fa",
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.promoter} -fo {output.promoter} -name
        bedtools getfasta -fi {input.genome} -bed {input.enhancer} -fo {output.enhancer} -name
        """   
