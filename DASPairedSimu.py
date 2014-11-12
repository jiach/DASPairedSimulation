#!/usr/bin/env python

import os
import sys
import simu_pairdas
from scipy.stats import uniform


# TODO: LOWER RANDOM NOISE
# TODO: REWRITE NULL PROPORTION -- FROM SAME VECTOR.

cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
pct_gene_picked = 20
num_sub = 40

# print set([x[0] for x in genes_picked[0]]).intersection(set([x[0] for x in genes_picked[1]]))
#print gene_arr.fetch_gene("IL1RN")

gene_anno_arr = simu_pairdas.RefSeqParser(cwd + "refgene_combined")
gene_arr = simu_pairdas.GeneArrReader(cwd + "allresult_filtered.txt")

# print simu_pairdas.fetch_gene_iso_comp(gene_arr.gene_arr, genes_picked[1])
#print "\n".join([str(x) for x in simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[1]), 10)])


#null_isocomp_arr = simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[0]), num_sub, 'null')

# print "\n".join([str(x[1][0]) for x in simu_pairdas.pick_null_genes(gene_arr.gene_arr)])

null_isocomp_arr = simu_pairdas.add_dither_to_isocomp(
    simu_pairdas.fetch_gene_iso_comp(gene_arr, simu_pairdas.pick_null_genes(gene_arr.gene_arr)), num_sub, 'null')



alt_isocomp_arr = simu_pairdas.add_dither_to_isocomp(
    simu_pairdas.fetch_gene_iso_comp(gene_arr, simu_pairdas.pick_alt_genes(gene_arr.gene_arr, pct_gene_picked)),
    num_sub, 'alt')
#print null_isocomp_arr['PDGFB']



null_fn = cwd+'simu_counts_null.txt'
null_fh = open(null_fn, 'w')
null_tot = len(null_isocomp_arr.keys())
index = 0
print("Output file NULL! Total: "+str(null_tot))
for gene_name in null_isocomp_arr.keys():
    sys.stdout.write(u"Progress: {0:d}   \r".format(index))
    sys.stdout.flush()
    index +=1
    for sub_idx in range(len(null_isocomp_arr[gene_name])):
        before = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), null_isocomp_arr[gene_name][sub_idx][1])
        after = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), null_isocomp_arr[gene_name][sub_idx][2])

        for i in range(len(before)):
            null_fh.write("\t".join([gene_name, str(sub_idx), str(before[i]), str(after[i])])+"\n")

null_fh.close()

alt_fn = cwd + 'simu_counts_alt.txt'
alt_fh = open(alt_fn, 'w')
alt_tot = len(alt_isocomp_arr.keys())
index = 0
print("Output file ALT! Total: "+str(alt_tot))
for gene_name in alt_isocomp_arr.keys():
    sys.stdout.write(u"Progress: {0:d}   \r".format(index))
    sys.stdout.flush()
    index +=1
    for sub_idx in range(len(alt_isocomp_arr[gene_name])):
        before = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000 + uniform.rvs(loc=-500, scale=1000, size=1),
                                                                    alt_isocomp_arr[gene_name][sub_idx][1])
        after = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000 + uniform.rvs(loc=-500, scale=1000, size=1),
                                                                   alt_isocomp_arr[gene_name][sub_idx][2])
        for i in range(len(before)):
            alt_fh.write("\t".join([gene_name, str(sub_idx), str(before[i]), str(after[i])]) + "\n")
alt_fh.close()