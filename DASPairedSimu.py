#!/usr/bin/env python

import simu_pairdas
import os

cwd = os.getenv("HOME") + "/Dropbox/Dissertation_2014/DAS_Paird/DASPairedSimulation/"
fname = cwd + "allresult_filtered.txt"
pct = 20

gene_arr = simu_pairdas.read_in_gene_arr(fname)
genes_picked = simu_pairdas.pick_null_alt_genes(gene_arr,pct)

for null_gene in genes_picked[0]:
    gene_arr[null_gene].