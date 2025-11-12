import os
include: "MetaInfo.smk"

rule mapping:
	input:
		expand("results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam", sample=SAMPLEs, rep=Reps),
		expand("results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",sample=SAMPLEs, rep=Reps),
		expand("results/uniq/{sample}_{rep}.uniq.bam", sample=SAMPLEs, rep=Reps),
		expand("results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.infer.txt", sample=SAMPLEs, rep=Reps),
		expand("results/bw_specific/{sample}_{rep}.Reverse.bw",sample=SAMPLEs, rep=Reps),

rule star_mapping:
	input:
		fastq_R1 = os.path.join(in_dir, "{sample}_{rep}_1.fq.gz"),
		fastq_R2 = os.path.join(in_dir, "{sample}_{rep}_2.fq.gz"),
	output:
		bam_R1 = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam",
	threads: 20
	params:
		index = STAR_INDX,
		gtf = Gh_GTF,
		outdir = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.",
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
		--alignIntronMin 1 \
		--scoreGapNoncan -4 \
		--scoreGapATAC -4 \
		--chimSegmentMin 15 \
		--chimJunctionOverhangMin 15 \
		--alignSJoverhangMin 15 \
		--alignSJDBoverhangMin 10 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--readFilesCommand gunzip -c 
		"""

rule sort:
	input:
		bam = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.bam",
	output:
		bam = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
	params:
		software = "samtools"
	threads: 10
	shell:"""
		{params.software} sort -@ {threads} {input.bam} -o {output.bam}
		{params.software} index -@ {threads} {output.bam}
		"""

rule uniq:
    input:
        bam =  "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
    output:
        bam = "results/uniq/{sample}_{rep}.uniq.bam",
        txt = "results/uniq/{sample}_{rep}_metrics.txt",
    log: "results/uniq/{sample}_{rep}.log",
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true  --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.txt} &> {log}
        """

rule infer_experiment:
    input:
        bam = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.Aligned.out.sort.bam",
        bed = Gh_BED,
    output:
        txt = "results/mapping/{sample}_{rep}_mapping/{sample}_{rep}.infer.txt"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}
        """

rule bam2wig_specific:
    input:
        bam = "results/uniq/{sample}_{rep}.uniq.bam",
        chromSize = CHROM_SIZE,
    output:
        bw = "results/bw_specific/{sample}_{rep}.Reverse.bw",
    threads:10
    params:
        prefix = "results/bw_specific/{sample}_{rep}",
    shell:
        """
		samtools index -@ {threads} {input.bam}
        bam2wig.py -i {input.bam} -s {input.chromSize} -d 1+-,1-+,2++,2-- -u -t 100000000 -o {params.prefix}
        """
