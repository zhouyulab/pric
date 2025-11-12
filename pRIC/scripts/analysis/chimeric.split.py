import os
import sys
from collections import defaultdict

GeneDict=defaultdict(int)
InterDict = defaultdict(int)
f_in = sys.argv[1]
f_out = sys.argv[2]
f_inter = sys.argv[3]

with open(f_in) as f:
    for line in f:
        if line.startswith("DonorChrom"):
            continue
        data = line.strip().split("\t")
        readnum, donorgene, acceptorgene, group = int(data[-4]), data[-3], data[-2], data[-1]
        if group != "intra-gene":
            continue
        donorgene = donorgene.split(";")
        acceptorgene = acceptorgene.split(";")
        gene_li = list(set(donorgene) & set(acceptorgene))
        for gene in gene_li:
            GeneDict[gene] += readnum

with open(f_in) as f:
    for line in f:
        if line.startswith("DonorChrom"):
            continue
        data = line.strip().split("\t")
        readnum, donorgene, acceptorgene, group = int(data[-4]), data[-3], data[-2], data[-1]
        if group != "inter-gene":
            continue
        donorgene = donorgene.split(";")
        acceptorgene = acceptorgene.split(";")
        for donor in donorgene:
            for accept in acceptorgene:
                if donor<accept:
                    key = (donor, accept)
                else:
                    key = (accept, donor)
                InterDict[key] += readnum


with open(f_out, "w") as f_o:
    f_o.write("GeneName\tReadNum\n")
    for gene in GeneDict:
        f_o.write(gene+"\t"+str(GeneDict[gene])+"\n")


with open(f_inter, "w") as f_o:
    f_o.write("Gene1\tGene2\tReadNum\n")
    for gene1, gene2 in InterDict:
        readnum = InterDict[(gene1, gene2)]
        f_o.write(gene1+"\t"+gene2+"\t"+str(readnum)+"\n")