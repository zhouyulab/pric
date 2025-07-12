import os
include: "MetaInfo_Tetraploid.smk"

rule rmDuprRNA:
	input:
		expand("../results/prepare/fastp/{sample}.R1.fq",sample=SAMPLEs),
		expand("../results/prepare/rmDuplicate/{sample}.rmDup.R1.fastq.gz",sample=SAMPLEs),
		expand("../results/prepare/fastqc_after/{sample}.R2_fastqc.html",sample=SAMPLEs),
		expand("../results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz", sample=SAMPLEs),

rule fastp:
    input:
        fq1 = os.path.join(data_indir, "{sample}_1.fq.gz"),
        fq2 = os.path.join(data_indir, "{sample}_2.fq.gz"),
    output:
        html = "results/prepare/fastp/{sample}.html",
        json = "results/prepare/fastp/{sample}.json",
        R1 = "results/prepare/fastp/{sample}.R1.fq",
        R2 = "results/prepare/fastp/{sample}.R2.fq",
    params:
        m = 20,
    threads: 8
    log: "results/prepare/fastp/{sample}.fastp.log",
    shell:
        """
        fastp -h {output.html} -j {output.json} -g -l {params.m} -w {threads} -i {input.fq1} -o {output.R1} -I {input.fq2} -O {output.R2} > {log} 2>&1
        """

rule fastqc:
	input:
		fq1 = "results/prepare/fastp/{sample}.R1.fq",
		fq2 = "results/prepare/fastp/{sample}.R2.fq",
	output:
		odir = "results/prepare/fastqc_after/{sample}.R2_fastqc.html",
	threads:
		8
	params:
		odir = "results/prepare/fastqc_after/",
	log:
		"results/prepare/fastqc/{sample}.log",
	shell:
		"""
		fastqc -o {params.odir} -t {threads} {input.fq1} &>{log}
        fastqc -o {params.odir} -t {threads} {input.fq2} &>{log}
		"""

rule trim:
    input:
        fq1 = os.path.join(data_indir, "{sample}_1.fq.gz"),
        fq2 = os.path.join(data_indir, "{sample}_2.fq.gz"),
    output:
        fq1 = "results/prepare/trim/{sample}/{sample}.R1_val_1.fq.gz",
        fq2 = "results/prepare/trim/{sample}/{sample}.R2_val_2.fq.gz",
    params:
        prefix = "results/prepare/trim/{sample}"
    threads: 8
    shell:
        """
        trim_galore --paired --quality 20 --length 15 -j {threads} --fastqc --small_rna -o {params.prefix} {input.fq1} {input.fq2}
        """

### remove_duplicated_reads.pl downloaded from https://github.com/caochch/RIC-seq.
rule rmDuplicate:
	input:
		fq_R1 = "results/prepare/fastp/{sample}.R1.fq",
		fq_R2 = "results/prepare/fastp/{sample}.R2.fq",
	output:
		R1 = "results/prepare/rmDuplicate/{sample}.rmDup.R1.fastq.gz",
		R2 = "results/prepare/rmDuplicate/{sample}.rmDup.R2.fastq.gz",
	params:
		remove_duplicated_reads = "perl remove_duplicated_reads.pl",
	threads: 8
	shell:
		"""
		{params.remove_duplicated_reads} {input.fq_R1} {input.fq_R2} {wildcards.sample}.R1.rmDup.fastq {wildcards.sample}.R2.rmDup.fastq 		
		gzip -c {wildcards.sample}.R1.rmDup.fastq > {output.R1}
		gzip -c {wildcards.sample}.R2.rmDup.fastq > {output.R2}
		rm {wildcards.sample}.R1.rmDup.fastq {wildcards.sample}.R2.rmDup.fastq
		"""

rule rm_rRNA:
	input:
		fastq_R1 = "results/prepare/rmDuplicate/{sample}.rmDup.R1.fastq.gz",
		fastq_R2 = "results/prepare/rmDuplicate/{sample}.rmDup.R2.fastq.gz",
	output:
		bam = "results/prepare/rmrRNA/{sample}_mapping/Aligned.sortedByCoord.out.bam",
		unmap_R1 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R1.fastq.gz",
		unmap_R2 = "results/prepare/rmrRNA/{sample}_mapping/{sample}.unmap.R2.fastq.gz",
	threads: 20
	params:
		index = rRNA_INDX,
		mapping_dir = "results/prepare/rmrRNA/{sample}_mapping",
	shell:"""
		if [[ -e {params.mapping_dir} ]]; then
			rm -r {params.mapping_dir}
		fi
		mkdir -p {params.mapping_dir}

		STAR --runThreadN {threads} \
		--genomeDir {params.index} \
		--outFileNamePrefix {params.mapping_dir}/ \
		--readFilesIn {input.fastq_R1} {input.fastq_R2} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 100 --outSAMattributes All \
		--alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
		--chimSegmentMin 15 --chimJunctionOverhangMin 15 \
		--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--outReadsUnmapped Fastx \
		--limitOutSJcollapsed 500000000 \
		--readFilesCommand gunzip -c --limitBAMsortRAM 23056060626
		
		gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
		gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
		rm {params.mapping_dir}/Unmapped.out.mate1 {params.mapping_dir}/Unmapped.out.mate2
		"""

