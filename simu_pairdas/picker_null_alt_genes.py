__author__ = 'cheng'
import operator
from scipy.stats import truncnorm
from parser_isoform_comp import *

# TODO: 1 gene repeated 1000 times

def pick_alt_genes(gene_arr, pct):
    """

   :rtype : list of lists of strings
   list[0] = ['name', ... ] for 0 - pct_picked_gene^th percentile of
the gene average hell d
   list[1] = ['name', ... ] for (100-pct_picked_gene)^th-100^th
percentile of the gene average hell d
   """
    gene_max_helld = {}

    for key, value in gene_arr.iteritems():
        min_max = value.getMinMaxHellingerDistance()
        gene_max_helld[key] = min_max['max']

    sorted_gene_by_max_helld = sorted(gene_max_helld.items(), key=operator.itemgetter(1))

    last_pct = sorted_gene_by_max_helld[-int(len(sorted_gene_by_max_helld) * pct / 100):]

    return last_pct


def pick_null_genes(gene_arr):
    null_picked = []
    for key, value in gene_arr.iteritems():
        null_picked.append((key, value.getRandomDistance()))

    return null_picked


def fetch_gene_iso_comp(gene_arr, genes_picked):
    """
   :rtype : dict
   :param gene_arr: gene_arr obj generated from read_in_gene_arr function,
                    a dict() with gene name as key, and
IsoformPropParser objects as value
   :param genes_picked: gene list[0] or list[1] generated from
pick_null_alt_genes function,
   :return: array of iso_comp
   """
    isocomp_arr = {}
    print "Fetch Gene Iso Comp!"
    for i in range(len(genes_picked)):
        isocomp = gene_arr.gene_arr[genes_picked[i][0]].getIsoComp(genes_picked[i][1][1])
        if isocomp is not None:
            isocomp_arr[genes_picked[i][0]] = isocomp
    return isocomp_arr


def add_dither_to_isocomp(isocomp_arr, num_subjects, type):
    """
    :param isocomp_arr:list of isoform compositions output from fetch_gene_iso_comp
    :param num_subjects: sample size that you need to multiply the original isocomp by.
    :return:  [[]]*num_subjects. inner [] =[[iso_name],[before],[after]]
    """
    tot_len = len(isocomp_arr.keys())
    print "Adding Dither to Isocomp! " + str(tot_len) + " Items!"
    dithered_isocomp_arr = {}
    index = 0
    for i in isocomp_arr.keys():
        sys.stdout.write(u"Progress: {0:d}   \r".format(index))
        sys.stdout.flush()
        dithered_isocomp_arr[i] = normalize_dither_compositions(isocomp_arr[i], num_subjects, type)
        index += 1
    return dithered_isocomp_arr


def normalize_dither_compositions(mat, num_subjects, type):
    arr_isocomp_bysub = [[]] * num_subjects

    for t in range(num_subjects):
        iso_name = [""] * len(mat[0])
        before = [0] * len(mat[0])
        after = [0] * len(mat[0])

        for i in range(len(mat[0])):
            iso_name[i] = mat[0][i]
            before[i] = mat[1][i] + truncnorm.rvs(loc=0, scale=0.8, a=-mat[1][i], b=1 - mat[1][i], size=1)[0]
            if type == "null":
                after[i] = mat[1][i] + truncnorm.rvs(loc=0, scale=0.8, a=-mat[1][i], b=1 - mat[1][i], size=1)[0]
            else:
                after[i] = mat[2][i] + truncnorm.rvs(loc=0, scale=0.8, a=-mat[2][i], b=1 - mat[2][i], size=1)[0]

        for i in range(len(mat[0])):
            before[i] = before[i] / sum(before)
            after[i] = after[i] / sum(after)
        arr_isocomp_bysub[t] = [iso_name, before, after]
    return arr_isocomp_bysub


class GeneArrReader:
    def __init__(self, fname):
        self.gene_arr = read_in_gene_arr(fname)

    def fetch_gene(self, gene_name):
        return self.gene_arr[gene_name]


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