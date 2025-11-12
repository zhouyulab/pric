import os
include: "MetaInfo.smk"

rule mapping:
	input:
		expand("./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam", sample=SAMPLEs, rep=Reps),
		expand("./results/uniq/{sample}_{rep}.uniq.bam", sample=SAMPLEs, rep=Reps),
		expand("./results/bw/{sample}_{rep}.bw",sample=SAMPLEs, rep=Reps),
		expand("./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.infer.txt",sample=SAMPLEs, rep=Reps),

rule star_mapping:
	input:
		fastq_R1 = "./results/prepare/rmrRNA/{sample}_{rep}_mapping/{sample}_{rep}.unmap.R1.fastq.gz",
		fastq_R2 = "./results/prepare/rmrRNA/{sample}_{rep}_mapping/{sample}_{rep}.unmap.R2.fastq.gz",
	output:
		bam = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam",
	threads: 20
	params:
		index = STAR_INDX,
		gtf = Gh_GTF,
		outdir = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.",
	shell:"""
		STAR --runMode alignReads \
		--runThreadN {threads} \
		--genomeDir {params.index} \
		--sjdbGTFfile {params.gtf} \
		--readFilesIn {input.fastq_R1} {input.fastq_R2} \
		--outFileNamePrefix {params.outdir} \
		--outSAMtype BAM Unsorted \
		--outReadsUnmapped Fastx \
		--outFilterMultimapNmax 100 \
		--outSAMattributes All \
		--readFilesCommand gunzip -c 
		"""

rule sort:
	input:
		bam = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam",
	output:
		bam = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
	params:
		software = "samtools"
	shell:"""
		{params.software} sort -@ 10 {input.bam} -o {output.bam}
		{params.software} index -@ 10 {output.bam}
		"""

rule uniq:
    input:
        bam =  "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
    output:
        bam = "./results/uniq/{sample}_{rep}.uniq.bam",
        txt = "./results/uniq/{sample}_{rep}_metrics.txt",
    log: "./results/uniq/{sample}_{rep}.log",
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true  --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.txt} &> {log}
        """

rule bam2wig:
    input:
        bam = "./results/uniq/{sample}_{rep}.uniq.bam",
        chromSize = CHROM_SIZE,
    output:
        bw = "./results/bw/{sample}_{rep}.bw",
    threads:10
    params:
        prefix = "./results/bw/{sample}_{rep}",
    shell:
        """
        samtools index -@ {threads} {input.bam}
        bam2wig.py -i {input.bam} -s {input.chromSize} -u -t 100000000 -o {params.prefix}
        """

rule infer_experiment:
    input:
        bam = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
        bed = Gh_BED,
    output:
        txt = "./results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.infer.txt"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}
        """
