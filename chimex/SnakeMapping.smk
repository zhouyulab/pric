include: "MetaInfo_Tetraploid.smk"

rule mapping:
	input:
		expand("results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam", sample=SAMPLEs),
		expand("results/mapping/{sample}/STAR_mapping/{sample}.Aligned.sort.bam",sample=SAMPLEs),
		expand("results/mapping/{sample}/STAR_mapping/{sample}.infer.txt",sample=SAMPLEs),

rule star_mapping:
	input:
		fastq_R1 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
		fastq_R2 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
	output:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam",
		sam = "results/mapping/{sample}/STAR_mapping/{sample}.Chimeric.out.sam",
	threads: 20
	params:
		index = STAR_INDX,
		gtf = Gh_GTF,
		outdir = "results/mapping/{sample}/STAR_mapping/{sample}.",
	shell:"""
		STAR --runMode alignReads \
		--runThreadN {threads} \
		--genomeDir {params.index} \
		--sjdbGTFfile {params.gtf} \
		--readFilesIn {input.fastq_R1} {input.fastq_R2} \
		--outFileNamePrefix {params.outdir} \
		--outSAMtype BAM Unsorted \
		--chimOutType SeparateSAMold \
		--outReadsUnmapped Fastx \
		--outFilterMultimapNmax 100 \
		--outSAMattributes All \
		--alignIntronMin 1 \
		--scoreGapNoncan -4 \
		--scoreGapATAC -4 \
		--chimSegmentMin 15 \
		--chimJunctionOverhangMin 15 \
		--outFilterMatchNminOverLread 0.33 \
		--outFilterScoreMinOverLread 0.33 \
		--alignSJoverhangMin 15 \
		--alignSJDBoverhangMin 10 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--readFilesCommand gunzip -c \
		--limitOutSJcollapsed 50000000
		"""

rule sort:
	input:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam",
	output:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.sort.bam",
	threads: 10
	shell:
		"""
		samtools sort -@ {threads} -o {output.bam} {input.bam}
		samtools index {output.bam}
		"""

rule infer_experiment:
	input:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam",
		Bed = Gh_BED,
	output:
		txt = "results/mapping/{sample}/STAR_mapping/{sample}.infer.txt",
	shell:
		"""
		infer_experiment.py -i {input.bam} -r {input.Bed} > {output.txt}
		"""

rule uniq:
	input:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam",
	output:
		bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.sort.uniq.bam",
	threads: 10
	shell:
		"""
		samtools view -@ {threads} -hb -q 30 -F 256 {input.bam} | samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		"""

rule not_fully_reads:
    input:
        fastq_R1 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
        fastq_R2 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
        bam = "results/mapping/{sample}/STAR_mapping/{sample}.Aligned.out.bam",
        sam = "results/mapping/{sample}/STAR_mapping/{sample}.Chimeric.out.sam",
    output:
        unmap_fq1 = "results/mapping/{sample}/BWA_mapping/{sample}_R1_toGenomeUnmapped_really.fq.gz",
        unmap_fq2 = "results/mapping/{sample}/BWA_mapping/{sample}_R2_toGenomeUnmapped_really.fq.gz",
    params:
        code_path = "python scripts/mapping/select_not_fully_aligned.py",
        unmap_fq1 = "results/mapping/{sample}/BWA_mapping/{sample}_R1_toGenomeUnmapped_really.fq",
        unmap_fq2 = "results/mapping/{sample}/BWA_mapping/{sample}_R2_toGenomeUnmapped_really.fq",
    shell:"""
        {params.code_path} -f1 {input.fastq_R1} -f2 {input.fastq_R2} -b {input.bam} -c {input.sam} -o1 {params.unmap_fq1} -o2 {params.unmap_fq2}
        gzip {params.unmap_fq1}
        gzip {params.unmap_fq2}
        """

# merge fastq with fastq-join
rule bwa_mapping:
	input:
		fq = "results/mapping/{sample}/BWA_mapping/{sample}_collapse.fq.gz",
	output:
		sam = "results/mapping/{sample}/BWA_mapping/{sample}.sam"
	threads:
		20
	params:
		index = BWA_INDX
	shell:
		"""
		bwa mem -t {threads} -k 12 -T 15 {params.index} {input.fq} > {output.sam}
		"""
