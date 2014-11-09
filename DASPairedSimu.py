#!/usr/bin/env python

import os
import simu_pairdas

cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
pct_gene_picked = 20
num_sub = 2
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
#x = 'PDGFB'
#print gene_anno_arr.gene_arr[x].getIsoReadCounts(1000000,null_isocomp_arr[x][1][1])

