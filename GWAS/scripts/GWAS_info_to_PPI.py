
Trait_dict = {
    'fiber maturity': 'fiber maturity',
    'FM': 'fiber maturity',
    "MAT": 'fiber maturity',
    'fiber length': 'fiber length',
    'fiber length uniformity': 'fiber length',
    'LU': 'fiber length', # fiber length uniformity 
    'FL': 'fiber length',
    'flowering date': 'flowering date',
    'FD': 'flowering date',
    'fiber strength': 'fiber strength',
    'FS': 'fiber strength',
    'LP': 'lint percentage',
    'lint index': 'lint percentage',
    'lint percentage': 'lint percentage',
    'plant height': 'plant height',
    'PH': 'plant height',
    'whole growth period': 'whole growth period',
    'LPA': 'leaf pubescense amount',
    'fiber micronaire': 'fiber micronaire',
    'M': 'fiber micronaire',
    'FE4': 'fiber elongation',
    'FE': 'fiber elongation',
    'E': 'fiber elongation',
    'fiber elongation': 'fiber elongation',
    'seed index': 'seed index',
    'boll weight': 'boll weight',
    'FWPB': "boll weight",
    'fiber quality': 'fiber quality',
    'spinning consistency index': 'other',
    'fiber weight per boll': 'other',
    'BN': 'other',
    'yield percentage before frost': 'other',
    'HNFFB': 'other', #'height of the node of the first fruiting branch', 
    'LI': "lint index",
    'BW': 'other',
    'SI': 'other',
    'SCI': 'other', #spinning consistency index 
    'height of the node of the first fruiting branch': 'other',
    
    }


import sys
from collections import defaultdict
f_in = sys.argv[1]
f_edge = sys.argv[2]
f_node = sys.argv[3]

trait_set = defaultdict(int)
with open(f_in) as f, open(f_edge, "w") as fe, open(f_node, "w") as fn:
    fe.write("\t".join(["Enhancer", "Promoter", "ReadNum", "Trait"])+"\n")
    fn.write("\t".join(["feature", "Chrom", "Group","GWAS","Trait", "GeneName", "AtGene", "GeneSymbol"])+"\n")
    for line in f:
        if line.startswith("Enhancer"):
            continue
        data = line.strip().split("\t")
        Enhancer,EnhancerChrom,Promoter,PromoterChrom,SNPloci, ReadNum, GeneName,AtGene,GeneSymbol,EnhancerSNP,PromoterSNP  = data
        EP_trait = set()
        if EnhancerSNP != "None":
            EnhancerGWAS = "GWAS"
            Enhancer_infos = EnhancerSNP.split("|")
            for Enhancer_info in Enhancer_infos:
                trait = Enhancer_info.split(":")[-1]
                EP_trait.add(trait)
        else:
            EnhancerGWAS = "None"


        if PromoterSNP != "None":
            PromoterGWAS = "GWAS"
            Promoter_infos = PromoterSNP.split("|")
            for Promoter_info in Promoter_infos:
                trait = Promoter_info.split(":")[-1]
                EP_trait.add(trait)
        else:
            PromoterGWAS = "None"
        EP_trait = list(set([Trait_dict[x] for x in EP_trait]))
        if len(list(EP_trait))>1:
            EP_trait = [x for x in EP_trait if x not in ['fiber quality', 'other']]
        EP_trait = EP_trait[0]
        trait_set[EP_trait] += 1
        edge_line = [Enhancer, Promoter,ReadNum, EP_trait]
        fe.write("\t".join(edge_line)+"\n")

        promoter_node_line = [Promoter, PromoterChrom, "Promoter", PromoterGWAS, EP_trait, GeneName, AtGene, GeneSymbol]
        enhancer_node_line = [Enhancer, EnhancerChrom, "Enhancer", EnhancerGWAS, EP_trait, Enhancer, Enhancer, Enhancer]
        fn.write("\t".join(promoter_node_line)+"\n")
        fn.write("\t".join(enhancer_node_line)+"\n")

