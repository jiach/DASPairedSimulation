#!/usr/bin/env python

import simu_pairdas
import os


cwd = os.getenv("HOME") + "/Dropbox/Dissertation_2014/DAS_Paird/DASPairedSimulation/"
fname = cwd + "allresult_filtered.txt"

gene_arr = {}
with open(fname, 'r') as fin:
    for line in fin:
        gene_name = line.split('\t')[0]
        if gene_name in gene_arr:
            gene_arr[gene_name].append(line)
        else:
            gene_arr[gene_name] = simu_pairdas.IsoformPropParser(line)



    # for key, value in gene_arr.iteritems():
    #     # hellDistList = value.getRenyiDiv(3)
    #     hellDistList = value.getHellingerDistance(False)
    #     for i in range(value.getNumSubjects()):
    #         print(key + "\t" + str(hellDistList[i]))

