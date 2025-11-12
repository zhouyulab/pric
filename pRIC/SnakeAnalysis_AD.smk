include: "MetaInfo_Tetraploid.smk"

merge_Sample = ["AD_WT_merge", "AD_fl_merge", "AD_pCp_merge"] 
bin_li = [5000000,2000000, 1000000, 500000, 200000, 100000, 10000, 50000, 1000, 2000]
ChIP_sample = {"AD_WT_merge": ["WT_H3K27ac_peaks.narrowPeak", "WT_H3K4me3_peaks.narrowPeak"],
               "AD_fl_merge": ["fl_H3K27ac_peaks.narrowPeak", "fl_H3K4me3_peaks.narrowPeak"],
               "AD_pCp_merge": ["H3K27ac_merge.peak", "H3K4me3_merge.peak"], }

rule RIC_chimeric:
    input:
        expand("./results/analysis/{sample}/{sample}.ChimericStatis.txt",sample = merge_Sample),
        expand("./results/analysis/{sample}/{sample}.intra.Chimeric.txt",sample = merge_Sample),
        expand("./results/analysis/{sample}/{sample}_{res}.bin.txt",sample = merge_Sample, res=bin_li),
        expand("./results/analysis/{sample}/{sample}.Intergenic.H3K27ac.txt",sample = merge_Sample),
        expand("./results/analysis/{sample}/{sample}.Intergenic.H3K4me3.txt",sample = merge_Sample),

rule group_chimeirc_reads:
    input:
        chimeric = "./results/chimeric/{sample}/{sample}.chimeric.cnt.txt",
        bed = "../Supplymental/Ghir_genome/merge_all.gene.bed6",
    output:
        txt =  "./results/analysis/{sample}/{sample}.ChimericGroup.txt",
        statis = "./results/analysis/{sample}/{sample}.ChimericStatis.txt",
    shell:
        """
        python scripts/analysis/chimeric.group.py -i {input.chimeric} -b {input.bed} -o {output.txt} -og {output.statis} 
        """

rule split_chimeric:
    input:
        txt =  "./results/analysis/{sample}/{sample}.ChimericGroup.txt",
    output:
        txt1 = "./results/analysis/{sample}/{sample}.intra.Chimeric.txt",
        txt2 = "./results/analysis/{sample}/{sample}.inter.Chimeric.txt",
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
        peak = lambda wildcards: "../AD_ChIP/results/peak/%s" % ChIP_sample[wildcards.sample][0],
        chimeric = "./results/analysis/{sample}/{sample}.ChimericGroup.txt",
    output:
        txt = "./results/analysis/{sample}/{sample}.Intergenic.H3K27ac.txt",
    shell:
        """
        python scripts/analysis/chimeric.H3K27ac.py -peak {input.peak} -c {input.chimeric} -o {output.txt}
        """

rule anno_intergenic_reads_in_peak_in_H3K4me3:
    input:
        peak = lambda wildcards: "../AD_ChIP/results/peak/%s" % ChIP_sample[wildcards.sample][1],
        chimeric = "./results/analysis/{sample}/{sample}.ChimericGroup.txt",
    output:
        txt = "./results/analysis/{sample}/{sample}.Intergenic.H3K4me3.txt",
    shell:
        """
        python scripts/analysis/chimeric.H3K4me4.py -peak {input.peak} -c {input.chimeric} -o {output.txt}
        """
