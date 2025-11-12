import sys
f_bed = sys.argv[1]
f_gtf = sys.argv[2]


with open(f_bed) as fi, open(f_gtf, "w") as fo:
    for line in fi:
        data = line.strip().split("\t")
        chrom, start, end, name, sample, strand, group = data
        start = str(int(start)+1)
        end = str(int(end)+1)
        gene_info = 'gene_type "%s"; gene_name "%s";'%("enhancer", name)
        trans_info = 'transcript_id "%s"; gene_type "%s"; gene_name "%s";'%(name, "enhancer", name)

        strand = "+"
        gtf_info_1 = [chrom, "enhancer" "gene", start, end, ".", strand, ".", gene_info]
        gtf_info_2 = [chrom, "enhancer", "transcript", start, end, ".", strand, ".", trans_info]
        gtf_info_3 = [chrom, "enhancer", "exon", start, end, ".", strand, ".", trans_info]
        for info in [gtf_info_1, gtf_info_2, gtf_info_3]:
            fo.write("\t".join(info) + "\n")
        strand = "-"
        gtf_info_1 = [chrom, "enhancer" "gene", start, end, ".", strand, ".", gene_info]
        gtf_info_2 = [chrom, "enhancer", "transcript", start, end, ".", strand, ".", trans_info]
        gtf_info_3 = [chrom, "enhancer", "exon", start, end, ".", strand, ".", trans_info]
        for info in [gtf_info_1, gtf_info_2, gtf_info_3]:
            fo.write("\t".join(info) + "\n")