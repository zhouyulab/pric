import os
include: "AD_MetaInfo.smk"

outdir = "./results/"
rule prepare:
    input:
        expand(os.path.join(outdir,"prepare/fastqc/{sample}_{ty}_{rep}_R2_fastqc.html"), sample=SAMPLEs, ty=Types, rep=Reps),

rule fastqc:
    input:
        fq1 = os.path.join(data_indir, "{sample}_{ty}_{rep}_R1.fq.gz"),
        fq2 = os.path.join(data_indir, "{sample}_{ty}_{rep}_R2.fq.gz"),
    output:
        os.path.join(outdir,"prepare/fastqc/{sample}_{ty}_{rep}_R2_fastqc.html")
    params:
        odir = os.path.join(outdir,"prepare/fastqc")
    log:
        os.path.join(outdir,"prepare/fastqc/{sample}_{ty}_{rep}.log")
    threads: 8
    shell:
        """
        fastqc -o {params.odir} -t {threads} {input.fq1} &>{log}
        fastqc -o {params.odir} -t {threads} {input.fq2} &>{log}
        """

