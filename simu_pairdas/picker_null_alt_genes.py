__author__ = 'cheng'
import operator

def pick_null_alt_genes(gene_arr, pct):
    """

    :rtype : list of lists of strings
    list[0] = ['name', ... ] for 0 - pct^th percentile of the gene average hell d
    list[1] = ['name', ... ] for (100-pct)^th-100^th percentile of the gene average hell d
    """
    gene_mean_helld = {}

    for key, value in gene_arr.iteritems():
            # hellDistList = value.getRenyiDiv(3)
            gene_mean_helld[key] = value.getAverageHellingerDistance()

    sorted_gene_by_mean_helld = sorted(gene_mean_helld.items(), key = operator.itemgetter(1))
    num_picked = int(len(sorted_gene_by_mean_helld)*pct/100)

    first_pct =[x[0] for x in sorted_gene_by_mean_helld[0:num_picked]]
    last_pct = [x[0] for x in sorted_gene_by_mean_helld[-num_picked:]]
    return([first_pct,last_pct])

def fetch_gene_iso_comp(gene_arr, genes_picked):
    pass