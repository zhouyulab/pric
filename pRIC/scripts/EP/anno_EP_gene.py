import sys
from collections import defaultdict

def record_anno(f_in):
    GeneAnno = defaultdict(list)
    with open(f_in) as fi:
        for line in fi:
            if line.startswith("cottonGene"):
                continue
            data = line.strip().split("\t")
            cotton_gene, ATgene, Symbol, Description, ifFiber, isC3 = data
            GeneAnno[cotton_gene].append([ATgene, Symbol, Description, ifFiber, isC3])
    return GeneAnno

f_EP = sys.argv[1]
f_anno = sys.argv[2]
f_out = sys.argv[3]

GeneAnno = record_anno(f_anno)
with open(f_EP) as fi, open(f_out, "w") as fo:
    for line in fi:
        data = line.strip().split("\t")
        if data[0] == "Enhancer":
            header = data + ["GeneName", "At", "Symbol", "Description", "isFiber", "isC3"]
            fo.write("\t".join(header)+"\n")
            continue
        promoter = data[2]
        gene_name = promoter.replace("pro", "")
        gene_info = GeneAnno[gene_name]
        if gene_info:
            outline = [*data, gene_name, *gene_info[0]]
        else:
            outline = [*data, gene_name] + ["None" for i in range(5)]
        fo.write("\t".join(outline)+"\n")
