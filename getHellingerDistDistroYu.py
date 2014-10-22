#!/usr/bin/env python

import IsoformPropParser
import os

cwd = os.getenv("HOME") + "/Dropbox/Dissertation_2014/DAS_Paird/DASPairedSimulation/"
fname = cwd + "allresult_filtered.txt"

geneArr = {}

with open(fname, 'r') as fin:
    for line in fin:
        geneName = line.split('\t')[0]
        if geneName in geneArr:
            geneArr[geneName].append(line)
        else:
            geneArr[geneName] = IsoformPropParser(line)

    for key, value in geneArr.iteritems():
        # hellDistList = value.getRenyiDiv(3)
        hellDistList = value.getHellingerDistance(False)
        for i in range(value.getNumSubjects()):
            print(key + "\t" + str(hellDistList[i]))
