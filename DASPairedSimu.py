#!/usr/bin/env python

import os
import simu_pairdas

cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
pct = 20

#print set([x[0] for x in genes_picked[0]]).intersection(set([x[0] for x in genes_picked[1]]))
#print gene_arr.fetch_gene("IL1RN")

# a = simu_pairdas.RefSeqParser(cwd + "refgene")
gene_arr = simu_pairdas.GeneArrReader(cwd + "allresult_filtered.txt")
genes_picked = simu_pairdas.pick_null_alt_genes(gene_arr.gene_arr, pct)
# print simu_pairdas.fetch_gene_iso_comp(gene_arr.gene_arr, genes_picked[1])
#print "\n".join([str(x) for x in simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[1]), 10)])

