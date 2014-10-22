__author__ = 'cheng'
import operator
from copy import  deepcopy

def pick_null_alt_genes(gene_arr, pct):
    """

    :rtype : list of lists of strings
    list[0] = ['name', ... ] for 0 - pct^th percentile of the gene average hell d
    list[1] = ['name', ... ] for (100-pct)^th-100^th percentile of the gene average hell d
    """
    gene_max_helld = {}
    gene_min_helld = {}

    for key, value in gene_arr.iteritems():
            # hellDistList = value.getRenyiDiv(3)
            min_max = value.getMinMaxHellingerDistance()
            gene_max_helld[key] = min_max['max']
            cur_min_helld = min_max['min']
            if cur_min_helld[0] >= 0:
                gene_min_helld[key] = cur_min_helld

    sorted_gene_by_max_helld = sorted(gene_max_helld.items(), key = operator.itemgetter(1))
    sorted_gene_by_min_helld = sorted(gene_min_helld.items(), key = operator.itemgetter(1))
    num_picked = int(len(sorted_gene_by_max_helld)*pct/100)

    first_pct =sorted_gene_by_min_helld[0:num_picked]
    last_pct = sorted_gene_by_max_helld[-num_picked:]
    return([first_pct,last_pct])

def fetch_gene_iso_comp(gene_arr, genes_picked):
    """
    :param gene_arr: gene_arr obj generated from read_in_gene_arr function,
                     a dict() with gene name as key, and IsoformPropParser objects as value
    :param genes_picked: gene list[0] or list[1] generated from pick_null_alt_genes function,
    :return: array of iso_comp
    """

    isocomp_arr = []
    for i in range(len(genes_picked)):
        isocomp_arr.append(gene_arr[genes_picked[i][0]].getIsoComp(genes_picked[i][1][1]))