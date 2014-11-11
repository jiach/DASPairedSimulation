#!/usr/bin/env python

import os
import simu_pairdas
from scipy.stats import uniform


# TODO: LOWER RANDOM NOISE
# TODO: REWRITE NULL PROPORTION -- FROM SAME VECTOR.

cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
pct_gene_picked = 20
num_sub = 40

#print set([x[0] for x in genes_picked[0]]).intersection(set([x[0] for x in genes_picked[1]]))
#print gene_arr.fetch_gene("IL1RN")

gene_anno_arr = simu_pairdas.RefSeqParser(cwd + "refgene_combined")
gene_arr = simu_pairdas.GeneArrReader(cwd + "allresult_filtered.txt")
genes_picked = simu_pairdas.pick_null_alt_genes(gene_arr.gene_arr, pct_gene_picked)

# print simu_pairdas.fetch_gene_iso_comp(gene_arr.gene_arr, genes_picked[1])
#print "\n".join([str(x) for x in simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[1]), 10)])


null_isocomp_arr = simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[0]), num_sub)
alt_isocomp_arr = simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[1]), num_sub)
#print null_isocomp_arr['PDGFB']

# for gene_name in null_isocomp_arr.keys():
#     for sub_idx in range(len(null_isocomp_arr[gene_name])):
#         before = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), null_isocomp_arr[gene_name][sub_idx][1])
#         after = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), null_isocomp_arr[gene_name][sub_idx][2])
#
#         for i in range(len(before)):
#             print "\t".join([gene_name, str(sub_idx), str(before[i]), str(after[i])])

for gene_name in alt_isocomp_arr.keys():
    for sub_idx in range(len(alt_isocomp_arr[gene_name])):
        before = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), alt_isocomp_arr[gene_name][sub_idx][1])
        after = gene_anno_arr.gene_arr[gene_name].getIsoReadCounts(5000+uniform.rvs(loc=-500, scale=1000, size=1), alt_isocomp_arr[gene_name][sub_idx][2])
        for i in range(len(before)):
            print "\t".join([gene_name, str(sub_idx), str(before[i]), str(after[i])])