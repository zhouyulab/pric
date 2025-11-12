import os
include: "MetaInfo.smk"

rule process_mapped_result:
	input:
		expand("./results/featureCount/{sample}_{rep}.protein.txt",  sample=SAMPLEs, rep=Reps),

rule featureCount:
	input:
		bam = "./results/uniq/{sample}_{rep}.uniq.bam",
		gtf = Gh_GTF,
	output:
		gene_txt = "./results/featureCount/{sample}_{rep}.protein.txt",
		summary = "./results/featureCount/{sample}_{rep}.protein.txt.summary"
	threads:
		10
	shell:
		"""
        featureCounts -s 2 -a {input.gtf} -o {output.gene_txt} -t exon -g gene_id --minOverlap 3 -p -C {input.bam} -T {threads}
		"""