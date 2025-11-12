import numpy as np
import os, sys

def record_percent(f_in):
    count_li = [0 for i in range(40)]
    total_num = 0
    with open(f_in) as f:
        for line in f:
            data = line.strip().split("\t")
            chimeric_count = data[4:]
            chimeric_count = list(map(int, chimeric_count))
            if len(chimeric_count) == 40:
                total_num += 1
                count_li = [count_li[i]+chimeric_count[i] for i in range(40)]
            else:
                print(line)
    precent_li = [count_li[i] * 100 /total_num for i in range(40)]
    return precent_li

f_eQTL = sys.argv[1]
f_random = sys.argv[2]
f_out = sys.argv[3]

eQTL_percent_li = record_percent(f_eQTL)
random_percent_li = record_percent(f_random)
header = ["Index", "Percent", "Group"]
with open(f_out, "w") as fo:
    fo.write("\t".join(header)+"\n")
    for indx, num in enumerate(eQTL_percent_li):
        outline = [indx, num, "eQTL"]
        outline = list(map(str, outline))
        fo.write("\t".join(outline)+"\n")
    for indx, num in enumerate(random_percent_li):
        outline = [indx, num, "random"]
        outline = list(map(str, outline))
        fo.write("\t".join(outline)+"\n")
