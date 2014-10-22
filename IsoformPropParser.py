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
             in range(len(mat))])
    else:
        return sum([1 / (1 - alpha) * (alpha * math.log(max(sys.float_info.min, mat[0][i]))) - (alpha - 1) * math.log(
            max(sys.float_info.min, mat[1][i])) for i in range(len(mat))])


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
        sumSqDiff = 0
        if len(isoformNames.difference(self.prop[0].keys())) == 0 and len(
                isoformNames.difference(self.prop[1].keys())) == 0:
            for isoform in list(isoformNames):
                sumSqDiff += (math.sqrt(self.prop[0][isoform]) - math.sqrt(self.prop[1][isoform])) ** 2
            return math.sqrt(sumSqDiff / 2)
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

    def getProp(self,isoformNames):
        return([[self.prop[0][x] for x in isoformNames],[self.prop[1][x] for x in isoformNames]])

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

    def getRenyiDiv(self, alpha):
        hellingerDist = []
        for value in self.subjectProp.values():
            hellD = value.getRenyiDiv(alpha, self.isoformList)
            hellingerDist.append(hellD)
        return hellingerDist

    def getAllProps(self):
        for key,value in 