__author__ = 'cheng'

__author__ = 'chengjia'

import sys
import math
import numpy as np


def calcRenyiDiv(mat, alpha):
    if alpha <= 0:
        return -1
    elif alpha == 1:
        return sum(
            [-2 * mat[0][i] * (math.log(max(sys.float_info.min, mat[0][i]) / max(sys.float_info.min, mat[1][i]))) for i
             in range(len(mat[0]))])
    else:
        return sum([1 / (1 - alpha) * (alpha * math.log(max(sys.float_info.min, mat[0][i]))) - (alpha - 1) * math.log(
            max(sys.float_info.min, mat[1][i])) for i in range(len(mat[0]))])

def calcHellingerDistance(mat):
    return math.sqrt(sum([(math.sqrt(mat[0][i])-math.sqrt(mat[1][i])) ** 2 for i in range(len(mat[0]))])/2)

def read_in_gene_arr(fname):
    gene_arr = {}
    with open(fname, 'r') as fin:
        for line in fin:
            gene_name = line.split('\t')[0]
            if gene_name in gene_arr:
                gene_arr[gene_name].append(line)
            else:
                gene_arr[gene_name] = IsoformPropParser(line)
    return gene_arr

class IsoformPropArr:
    def __init__(self, isoformName, groupIdx, prop):
        self.prop = {0: {}, 1: {}}
        self.prop[groupIdx - 1][isoformName] = prop

    def append(self, isoformName, groupIdx, prop):
        self.prop[groupIdx - 1][isoformName] = prop

    def verifySelf(self, isoformNames):
        if len(isoformNames.difference(self.prop[0].keys())) == 0 and len(
                isoformNames.difference(self.prop[1].keys())) == 0:
            if abs(sum(self.prop[0].values()) - 1) <= np.spacing(1) and abs(
                            sum(self.prop[1].values()) - 1) <= np.spacing(1):
                return True
            else:
                return False
        else:
            return False

    def getHellingerDistance(self, isoformNames):
        sum_before = sum(self.prop[0].values())
        sum_after = sum(self.prop[1].values())
        prop_mat = [[], []]
        if len(isoformNames.difference(self.prop[0].keys())) == 0 and len(
                isoformNames.difference(self.prop[1].keys())) == 0 and sum_before > 0.5 and sum_after > 0.5:
            for isoform in list(isoformNames):
                prop_mat[0].append(self.prop[0][isoform])
                prop_mat[1].append(self.prop[1][isoform])
            return calcHellingerDistance(prop_mat)
        else:
            return -1

    def printAll(self):
        print(self.prop)

    def getRenyiDiv(self, alpha, isoformNames):
        if len(isoformNames.difference(self.prop[0].keys())) == 0 and len(
                isoformNames.difference(self.prop[1].keys())) == 0:
            mat = [[], []]
            for isoform in list(isoformNames):
                mat[0].append(self.prop[0][isoform])
                mat[1].append(self.prop[1][isoform])
            return calcRenyiDiv(mat, alpha)
        else:
            return -1

    def getProp(self, isoformNames):
        iso_prop_mat = ([], [], [])
        isoformNames = list(isoformNames)
        for i in range(len(isoformNames)):
            iso_prop_mat[0].append(isoformNames[i])
            iso_prop_mat[1].append(self.prop[0][isoformNames[i]])
            iso_prop_mat[2].append(self.prop[1][isoformNames[i]])
        return iso_prop_mat

class IsoformPropParser:
    def __init__(self, line):
        lineFields = line.split('\t')
        self.geneName = lineFields[0]
        self.geneLocation = lineFields[1]
        self.isoformList = set()
        self.isoformList.add(lineFields[2])

        if int(lineFields[4]) > 7:
            subjectIdx = int(lineFields[4]) - 8
        else:
            subjectIdx = int(lineFields[4])

        self.subjectProp = {subjectIdx: IsoformPropArr(lineFields[2], int(lineFields[6]), float(lineFields[5]))}

    def append(self, line):

        lineFields = line.split('\t')
        self.isoformList.add(lineFields[2])

        if int(lineFields[4]) > 7:
            subjectIdx = int(lineFields[4]) - 8
        else:
            subjectIdx = int(lineFields[4])

        if subjectIdx in self.subjectProp:
            self.subjectProp[subjectIdx].append(lineFields[2], int(lineFields[6]), float(lineFields[5]))
        else:
            self.subjectProp[subjectIdx] = IsoformPropArr(lineFields[2], int(lineFields[6]), float(lineFields[5]))

    def getNumSubjects(self):
        return len(self.subjectProp)

    def printAll(self):
        for key, value in self.subjectProp.iteritems():
            print(str(key) + ":")
            value.printAll()

    def verifyAll(self):
        for key, value in self.subjectProp.iteritems():
            print(str(key) + ":")
            print(value.verifySelf(self.isoformList))

    def getHellingerDistance(self, verbose):
        hellingerDist = []
        for value in self.subjectProp.values():
            hellD = value.getHellingerDistance(self.isoformList)
            hellingerDist.append(hellD)
            if verbose:
                value.printAll()
                print(hellD)
        return hellingerDist

    def getAverageHellingerDistance(self):
        return np.mean(self.getHellingerDistance(False))

    def getMinMaxHellingerDistance(self):
        all_hell_darr = self.getHellingerDistance(False)
        non_neg_hell_darr = []
        non_neg = False
        min_max = dict()
        for x in all_hell_darr:
            if x >= 0:
                non_neg_hell_darr.append(x)
                non_neg = True
        if non_neg:
            min_non_neg_hell_d = np.min(non_neg_hell_darr)
            min_max['min'] = (min_non_neg_hell_d, all_hell_darr.index(min_non_neg_hell_d))
        else:
            min_max['min'] = (-1, -1)

        min_max['max'] = (np.max(all_hell_darr), np.argmax(all_hell_darr))

        return min_max

    def getRenyiDiv(self, alpha):
        hellingerDist = []
        for value in self.subjectProp.values():
            hellD = value.getRenyiDiv(alpha, self.isoformList)
            hellingerDist.append(hellD)
        return hellingerDist

    def getIsoComp(self, subjectIdx):
        if subjectIdx > self.getNumSubjects():
            return None
        else:
            return self.subjectProp[subjectIdx].getProp(self.isoformList)