include: "AD_MetaInfo.smk"

rule AD_feature: 
    input:
        "./results/peak/H3K27ac_merge.peak",
        "./results/peak/H3K4me3_merge.peak",
        "./results/peak/H3K4me3_Promoter.peak",
        "./results/feature/enhancer.bed",
        "./results/feature/promoter.bed",
        "./results/feature/enhancer.count.bed",
        "./results/feature/enhancer.ROSE.pdf",
        "./results/feature/enhancer.reads.norm.txt",
        "./results/feature/enhancer_nearest_promoter.txt",
        "./results/feature/Enhancer.region.pdf",
        "./results/feature/Promoter.region.pdf",
        "./results/feature/enhancer.group.anno.bed",
        "./results/feature/enhancer.group.fa",

rule merge_H3K27ac_peak:
    input:
        WT = "./results/peak/WT_H3K27ac_peaks.narrowPeak",
        fl = "./results/peak/fl_H3K27ac_peaks.narrowPeak",
    output:
        WT = "./results/peak/WT_H3K27ac_merge.peak",
        fl = "./results/peak/fl_H3K27ac_merge.peak",
        merge = "./results/peak/H3K27ac_merge.peak",
    shell:
        """
        python scripts/peaks/merge_peak.py --input-WT {input.WT} --input-fl {input.fl} --output-WT {output.WT} --output-fl {output.fl} -o {output.merge}
        """

rule merge_H3K4me3_peak:
    input:
        WT = "./results/peak/WT_H3K4me3_peaks.narrowPeak",
        fl = "./results/peak/fl_H3K4me3_peaks.narrowPeak",
    output:
        WT = "./results/peak/WT_H3K4me3_merge.peak",
        fl = "./results/peak/fl_H3K4me3_merge.peak",
        merge = "./results/peak/H3K4me3_merge.peak",
    shell:
        """
        python scripts/peaks/merge_peak.py --input-WT {input.WT} --input-fl {input.fl} --output-WT {output.WT} --output-fl {output.fl} -o {output.merge}
        """

rule anno_H3K4me3_peak:
    input:
        peak = "./results/peak/H3K4me3_merge.peak",
        Bed6 = GENE,
    output:
        txt = "./results/peak/H3K4me3_Promoter.peak",
    params:
        flank = 2000
    shell:
        """
        python scripts/peaks/anno_H3K4me3_peak.py -p {input.peak} -b {input.Bed6} -o {output.txt} -f {params.flank}
        """


rule merge_H3K27ac_with_ATAC_DNase_GRO:
    input:
        H3K27ac = "./results/peak/H3K27ac_merge.peak",
        DNase = "../DNase/results/peak/DNase_merge.peak",
        H3K4me4 = "./results/peak/H3K4me3_Promoter.peak",
        ATAC = "../scATAC/results/peak/ATAC_merge.peak",
        GRO = "../GRO/results/peak/AD_GRO.filter.bed",
    params:
        prefix = "./results/feature/"
    output:
        enhancer = "./results/feature/enhancer.bed",
        promoter = "./results/feature/promoter.bed", 
    shell:
        """
        python scripts/peaks/merge_peak_with_ATAC_DNase_GRO.py --DNase {input.DNase} --ATAC {input.ATAC} --H3K4me4 {input.H3K4me4} \
        --H3K27ac {input.H3K27ac} --GRO {input.GRO} -e {output.enhancer} -p {output.promoter} --prefix {params.prefix}
        """

rule enhancer_sigle:
    input:
        peak = "./results/feature/enhancer.bed",
        bam_rep1 = "./results/uniq/WT_H3K27ac_merge.bam",
        bam_rep2 = "./results/uniq/fl_H3K27ac_merge.bam",
        DNase_rep1 = "../DNase/results/uniq/DNase_Rep1.uniq.bam",
        DNase_rep2 = "../DNase/results/uniq/DNase_Rep2.uniq.bam",
    output:
        bed = "./results/feature/enhancer.count.bed"
    shell:
        """
        bedtools multicov -bed {input.peak} -bams {input.bam_rep1} {input.bam_rep2} {input.DNase_rep1} {input.DNase_rep2} > {output.bed}
        """

rule class_enhancer:
    input:
        bed = "./results/feature/enhancer.count.bed"
    output:
        bed = "./results/feature/enhancer.group.bed",
        pdf = "./results/feature/enhancer.ROSE.pdf",
    shell:
        """
        Rscript scripts/peaks/class_enhancr_to_TE_and_SE.R {input.bed} {output.bed} {output.pdf}
        """

rule enhancer_read_norm:
    input:
        enhancer = "./results/feature/enhancer.group.bed",
        RNA = expand("../RNA-seq/results/uniq/{sample}_{rep}.uniq.bam", sample=["WT", "fl"], rep=["Rep1", "Rep2"]),
        GRO = expand("../GRO/results/uniq/{samp}.uniq.bam", samp=["GRO_AD_WT_Rep1", "GRO_AD_WT_Rep2"]),
    output:
        gtf = "./results/feature/enhancer.group.gtf",
        txt = "./results/feature/enhancer.conunt.txt",
        summary = "./results/feature/enhancer.conunt.txt.summary",
        reads = "./results/feature/enhancer.reads.norm.txt",
    threads: 8
    shell:
        """
        python scripts/peaks/enhancer_bed_to_gtf.py {input.enhancer} {output.gtf}
        featureCounts -a {output.gtf} -o {output.txt}  -t exon -g gene_name --minOverlap 3 -O -p -C -T {threads} {input.RNA} {input.GRO}
        python scripts/peaks/enhancer_norm_reads.py {output.txt} {output.summary} {output.reads}
        """

rule find_enhancer_nearst_promoter:
    input:
        promoter = "./results/feature/promoter.bed",
        enhancer = "./results/feature/enhancer.group.bed",
    output:
        txt = "./results/feature/enhancer_nearest_promoter.txt",
    shell:
        """
        python scripts/peaks/find_enhancer_nearst_promoter.py -p {input.promoter} -e {input.enhancer} -o {output.txt}
        """

rule enhancer_feature:
    input:
        enhancer = "./results/feature/enhancer.group.bed",
        H3K27ac = "./results/bw/WT_H3K27ac.bw",
        H3K4me4 = "./results/bw/WT_H3K4me3.bw",
        DNase = "../DNase/results/bw/DNase_Rep1.bw",
        GRO_forward = "../GRO/results/bw/GRO_AD_WT_Rep2.Forward.bw",
        GRO_reverse = "../GRO/results/bw/GRO_AD_WT_Rep2.Reverse.bw",
        ATAC = "../scATAC/results/bw/WT.bw"
    output:
        txt = "./results/feature/Enhancer.region.gz",
        matrix = "./results/feature/Enhancer.outFileNameMatrix",
    threads: 8
    params:
        flank = 1000
    shell:
        """
        computeMatrix scale-regions --regionsFileName {input.enhancer} --scoreFileName {input.H3K27ac} {input.H3K4me4} {input.DNase} {input.GRO_forward} {input.GRO_reverse} {input.ATAC} \
        --outFileNameMatrix {output.matrix} --beforeRegionStartLength {params.flank} --afterRegionStartLength  {params.flank} --regionBodyLength 1000 --missingDataAsZero \
        --startLabel TSS --endLabel TES  --skipZeros --numberOfProcessors {threads} --outFileName {output.txt}
        """

rule plot_enhancer:
    input:
        txt = "./results/feature/Enhancer.region.gz",
    output:
        png = "./results/feature/Enhancer.region.pdf",
        txt = "./results/feature/Enhancer.region.Matrix.gz",
    shell:
        """
        plotHeatmap -m {input.txt} -out {output.png} --colorMap Blues Greens Oranges RdBu RdBu Reds \
        --zMin 0 0 0 -0.2 -0.2 0 --zMax 0.8 0.8 0.2 0.2 0.2 0.1 --heatmapHeight 15 --heatmapWidth 6 --outFileNameMatrix {output.txt}
        """

rule promoter_feature:
    input:
        promoter = "./results/feature/promoter.bed",
        H3K27ac = "./results/bw/WT_H3K27ac.bw",
        H3K4me4 = "./results/bw/WT_H3K4me3.bw",
        DNase = "../DNase/results/bw/DNase_Rep1.bw",
        GRO_forward = "../GRO/results/bw/GRO_AD_WT_Rep2.Forward.bw",
        GRO_reverse = "../GRO/results/bw/GRO_AD_WT_Rep2.Reverse.bw",
        ATAC = "../scATAC/results/bw/WT.bw"
    output:
        txt = "./results/feature/Promoter.region.gz",
        matrix = "./results/feature/Promoter.outFileNameMatrix",
    threads: 8
    params:
        flank = 1000
    shell:
        """
        computeMatrix scale-regions --regionsFileName {input.promoter} --scoreFileName {input.H3K27ac} {input.H3K4me4} {input.DNase} {input.GRO_forward} {input.GRO_reverse} {input.ATAC} \
        --outFileNameMatrix {output.matrix} --beforeRegionStartLength {params.flank} --afterRegionStartLength  {params.flank} --regionBodyLength 1000 --missingDataAsZero \
        --startLabel TSS --endLabel TES  --skipZeros --numberOfProcessors {threads} --outFileName {output.txt}
        """

rule plot_promoter:
    input:
        txt = "./results/feature/Promoter.region.gz",
    output:
        png = "./results/feature/Promoter.region.pdf",
        txt = "./results/feature/Promoter.region.Matrix.gz",
    shell:
        """
        plotHeatmap -m {input.txt} -out {output.png} --colorMap Blues Greens Oranges RdBu RdBu Reds \
        --zMin 0 0 0 -0.2 -0.2 0 --zMax 1.5 1.5 0.1 0.2 0.2 0.1 --heatmapHeight 15 --heatmapWidth 6 --outFileNameMatrix {output.txt}
        """

rule anno_enhancer:
    input:
        enhancer = "./results/feature/enhancer.group.bed",
        bed = GENE,
    output:
        txt = "./results/feature/enhancer.group.anno.bed",
    shell:
        """
        python scripts/peaks/anno_peak.py -r {input.bed} -e {input.enhancer} -o {output.txt}
        """

rule feature_fasta:
    input:
        promoter = "./results/feature/promoter.bed",
        enhancer = "./results/feature/enhancer.group.bed",
        genome = fasta,
    output:
        promoter = "./results/feature/promoter.fa",
        enhancer = "./results/feature/enhancer.group.fa",
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.promoter} -fo {output.promoter} -name
        bedtools getfasta -fi {input.genome} -bed {input.enhancer} -fo {output.enhancer} -name
        """   