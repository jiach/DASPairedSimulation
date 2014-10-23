#!/usr/bin/env python

import simu_pairdas
import os


cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
fname = cwd + "allresult_filtered.txt"
pct = 20

gene_arr = simu_pairdas.geneArrReader(fname)
genes_picked = simu_pairdas.pick_null_alt_genes(gene_arr.gene_arr, pct)


#print(simu_pairdas.fetch_gene_iso_comp(gene_arr, genes_picked[1]))
#print "\n".join([str(simu_pairdas.calcHellingerDistance([x[1], x[2]])) for x in simu_pairdas.add_dither_to_isocomp(simu_pairdas.fetch_gene_iso_comp(gene_arr,genes_picked[0]))])

#print set([x[0] for x in genes_picked[0]]).intersection(set([x[0] for x in genes_picked[1]]))
print gene_arr.fetch_gene("IL1RN")