import os
include: "MetaInfo.smk"

out_dir = "./results"
rule prepare:
	input:
		expand(os.path.join(out_dir, "prepare/fastqc","{sample}_{rep}_1_fastqc.html"), sample=SAMPLEs, rep=Reps),
		expand(os.path.join(out_dir, "prepare/rmrRNA/{sample}_{rep}_mapping/{sample}_{rep}.unmap.R1.fastq.gz"),sample=SAMPLEs, rep=Reps),

rule fastqc:
	input:
		fq_R1 = os.path.join(in_dir, "{sample}_{rep}_1.fq.gz"),
		fq_R2 = os.path.join(in_dir, "{sample}_{rep}_2.fq.gz"),
	output:
		fq_R1 = os.path.join(out_dir, "prepare/fastqc/","{sample}_{rep}_1_fastqc.html"),
		fq_R2 = os.path.join(out_dir, "prepare/fastqc/","{sample}_{rep}_2_fastqc.html"),
	threads:5
	log: os.path.join(out_dir, "prepare/fastqc/", "{sample}_{rep}.log"),
	params:
		outdir = os.path.join(out_dir, "prepare/fastqc")
	shell:
		"""
		fastqc -t {threads} -o {params.outdir} {input.fq_R1} &> {log}
		fastqc -t {threads} -o {params.outdir} {input.fq_R2} &> {log}
		"""


rule rm_rRNA:
	input:
		fq_R1 = os.path.join(in_dir, "{sample}_{rep}_1.fq.gz"),
		fq_R2 = os.path.join(in_dir, "{sample}_{rep}_2.fq.gz"),
	output:
		bam = os.path.join(out_dir, "prepare/rmrRNA/{sample}_{rep}_mapping/Aligned.sortedByCoord.out.bam"),
		unmap_R1 = os.path.join(out_dir, "prepare/rmrRNA/{sample}_{rep}_mapping/{sample}_{rep}.unmap.R1.fastq.gz"),
		unmap_R2 = os.path.join(out_dir, "prepare/rmrRNA/{sample}_{rep}_mapping/{sample}_{rep}.unmap.R2.fastq.gz"),
	threads: 24
	params:
		index = rRNA_INDX,
		mapping_dir = os.path.join(out_dir, "prepare/rmrRNA/{sample}_{rep}_mapping")
	shell:"""
		if [[ -e {params.mapping_dir} ]]; then
			rm -r {params.mapping_dir}
		fi
		mkdir -p {params.mapping_dir}

		STAR --runThreadN {threads} \
		--genomeDir {params.index} \
		--outFileNamePrefix {params.mapping_dir}/ \
		--readFilesIn {input.fq_R1} {input.fq_R2} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 100 --outSAMattributes All \
		--alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 \
		--chimSegmentMin 15 --chimJunctionOverhangMin 15 \
		--alignSJoverhangMin 15 --alignSJDBoverhangMin 10 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--outReadsUnmapped Fastx \
		--limitOutSJcollapsed 500000000 \
		--limitBAMsortRAM 23056060626 \
		--readFilesCommand gunzip -c 
		
		gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
		gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
		rm {params.mapping_dir}/Unmapped.out.mate1 {params.mapping_dir}/Unmapped.out.mate2
		"""
