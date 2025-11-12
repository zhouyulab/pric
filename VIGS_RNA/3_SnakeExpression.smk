import os
include: "MetaInfo.smk"

rule process_mapped_result:
	input:
		expand("results/featureCount/{sample}_{rep}.txt",sample=SAMPLEs, rep=Reps),
		expand("results/featureCount/{sample}_{rep}.FPKM.txt",sample=SAMPLEs, rep=Reps),
		expand("results/lncRNA/{sample}_{rep}_lncRNA.FPKM.txt",sample=SAMPLEs, rep=Reps),
		expand("results/merge/{sample}_{rep}.txt",sample=SAMPLEs, rep=Reps),

rule featureCount:
	input:
		bam = "results/uniq/{sample}_{rep}.uniq.bam",
		gtf = Gh_GTF,
	output:
		gene_txt = "results/featureCount/{sample}_{rep}.txt",
		summary = "results/featureCount/{sample}_{rep}.txt.summary",
	threads:
		4
	shell:
		"""
        featureCounts -s 2 -a {input.gtf} -o {output.gene_txt} -t exon -g gene_id --minOverlap 3 -p -C {input.bam} -T {threads}
		"""

rule calculte_FPKM:
	input:
		gene_txt = "results/featureCount/{sample}_{rep}.txt",
		summary = "results/featureCount/{sample}_{rep}.txt.summary",
	output:
		txt = "results/featureCount/{sample}_{rep}.FPKM.txt",
	shell:
		"""
		python scripts/calculate_gene_fpkm.py -i {input.gene_txt} -s {input.summary} -o {output.txt}
		"""

rule featureCount_lnCRNA:
	input:
		bam = "results/uniq/{sample}_{rep}.uniq.bam",
		gtf = "../Supplymental/genome/Ghir_genome/lnc.gtf",
	output:
		gene_txt = "results/lncRNA/{sample}_{rep}_lncRNA.txt",
		summary = "results/lncRNA/{sample}_{rep}_lncRNA.txt.summary",
	threads:
		4
	shell:
		"""
        featureCounts -s 2 -a {input.gtf} -o {output.gene_txt} -t exon -g gene_id --minOverlap 3 -p -C {input.bam} -T {threads}
		"""

rule calculte_FPKM_lncRNA:
	input:
		gene_txt = "results/lncRNA/{sample}_{rep}_lncRNA.txt",
		summary = "results/lncRNA/{sample}_{rep}_lncRNA.txt.summary",
	output:
		txt = "results/lncRNA/{sample}_{rep}_lncRNA.FPKM.txt",
	shell:
		"""
		python scripts/calculate_gene_fpkm.py -i {input.gene_txt} -s {input.summary} -o {output.txt}
		"""
