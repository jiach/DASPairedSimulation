__author__ = 'cheng'
import numpy as np

class GeneAnno:
    def __init__(self, name, contig, strand, start, end, isoforms):
        self.iso_len = []
        self.name = name
        self.contig = contig
        self.strand = strand
        self.start = start
        self.end = end
        self.isoforms = isoforms
        self.len = abs(self.start - self.end) + 1
        self.num_isoforms = len(self.isoforms)
        self.comp_mat = []
        self.region_len = []
        self.region_start = []
        self.region_end = []
        self.iso_major_comp_mat = []

    def addRegion(self, start, end, presence_isoform):
        self.region_start.append(start)
        self.region_end.append(end)
        self.region_len.append(abs(start - end) + 1)
        self.comp_mat.append(presence_isoform)

    def __str__(self):
        return "\n".join(["\t".join(
            [self.name, self.contig, self.strand, str(self.region_start[x]), str(self.region_end[x]),
             str(self.region_len[x]), ",".join([str(y) for y in self.comp_mat[x]])]) for x in range(len(self.region_len))])

    def getIsoMajorCompMat(self):
        self.iso_major_comp_mat = np.asarray(self.comp_mat).T

    def getIsoformLen(self):
        self.getIsoMajorCompMat()
        self.iso_len = np.dot(self.iso_major_comp_mat, np.asarray(self.region_len))

    def getIsoReadCounts(self, tot_reads, iso_prop):
        self.getIsoformLen()


class RefSeqParser:
    def __init__(self, fname):
        self.gene_arr = {}
        with open(fname, 'r') as fh:
            for line in fh.readlines():
                line_elems = line.split("\t")
                if not self.gene_arr.has_key(line_elems[0]):
                    self.gene_arr[line_elems[0]] = GeneAnno(line_elems[0], line_elems[1], line_elems[2],
                                                            int(line_elems[3]), int(line_elems[4]),
                                                            line_elems[5].rstrip(",\n").split(","))
                else:
                    self.gene_arr[line_elems[0]].addRegion(int(line_elems[3]), int(line_elems[4]),
                                                           [int(x) for x in line_elems[5].rstrip(",\n").split(",")])
        fh.close()

    def __str__(self):
        return "\n".join([str(x) for x in self.gene_arr.values()])