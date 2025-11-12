include: "MetaInfo_Diploid.smk"
res_li = [5000000,2000000, 1000000, 500000, 200000, 100000, 10000, 50000, 1000, 2000]

rule RIC_chimeric:
    input:
        expand("./results/analysis/{sample}.ChimericStatis.txt", sample=SAMPLEs+["A_WT_merge"]),
        expand("./results/analysis/{sample}.intra.Chimeric.txt", sample=SAMPLEs+["A_WT_merge"]),
        expand("./results/analysis/{sample}/{sample}_{res}.bin.txt",sample = SAMPLEs+["A_WT_merge"], res=res_li),
        expand("./results/analysis/{sample}.Intergenic.H3K27ac.txt", sample=SAMPLEs+["A_WT_merge"]),
        expand("./results/analysis/{sample}.Intergenic.H3K4me3.txt", sample=SAMPLEs+["A_WT_merge"]),

rule group_chimeirc_reads:
    input:
        chimeric = "./results/chimeric/{sample}/{sample}.chimeric.cnt.txt",
        bed = GENE,
    output:
        txt =  "./results/analysis/{sample}.ChimericGroup.txt",
        statis = "./results/analysis/{sample}.ChimericStatis.txt",
    shell:
        """
        python scripts/analysis/chimeric.group.py -i {input.chimeric} -b {input.bed} -o {output.txt} -og {output.statis} 
        """

rule split_chimeric:
    input:
        txt = "./results/analysis/{sample}.ChimericGroup.txt",
    output:
        txt1 = "./results/analysis/{sample}.intra.Chimeric.txt",
        txt2 = "./results/analysis/{sample}.inter.Chimeric.txt",
    shell:
        """
        python scripts/analysis/chimeric.split.py {input.txt} {output.txt1} {output.txt2}
        """

rule global_chimeric_bin:
    input:
        chimeric = "./results/chimeric/{sample}/{sample}.chimeric.cnt.txt",
        chrom_size = CHROM_SIZE,
    params:
        resolution = "{res}"
    output:
        txt = "./results/analysis/{sample}/{sample}_{res}.bin.txt",
        id = "./results/analysis/{sample}/{sample}_{res}.ID.txt",
    shell:
        """
        python scripts/analysis/chimeric.bin.py -i {input.chimeric} -c {input.chrom_size} -r {params.resolution} -o {output.txt} --out-id {output.id}
        """

rule anno_intergenic_reads_in_peak_in_H3K27ac:
    input:
        peak = "../AA_ChIP/results/peak/H3K4me3_peaks.narrowPeak",
        chimeric = "./results/analysis/{sample}.ChimericGroup.txt",
    output:
        txt = "./results/analysis/{sample}.Intergenic.H3K27ac.txt",
    shell:
        """
        python scripts/analysis/chimeric.H3K27ac.py -peak {input.peak} -c {input.chimeric} -o {output.txt}
        """

rule anno_intergenic_reads_in_peak_in_H3K4me3:
    input:
        peak = "../AA_ChIP/results/peak/H3K27ac_peaks.narrowPeak",
        chimeric = "./results/analysis/{sample}.ChimericGroup.txt",
    output:
        txt = "./results/analysis/{sample}.Intergenic.H3K4me3.txt",
    shell:
        """
        python scripts/analysis/chimeric.H3K4me4.py -peak {input.peak} -c {input.chimeric} -o {output.txt}
        """
