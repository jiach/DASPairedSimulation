#!/usr/bin/env python
import os

__author__ = 'chengjia'

cwd = os.getenv("HOME") + "/IdeaProjects/DASPairedSimulation/"
fn = cwd + 'refgene_ori'

fh = open(fn, 'r')
lines = fh.readlines()

print "".join(["\t".join([x.split("\t")[0]+"_duplicate"]+x.split("\t")[1:6]) for x in lines])
#for line in lines:
    # if not line.split("\t")[0].endswith("duplicate"):
    #     print line.rstrip("\n")