__author__ = 'cheng'


class GeneAnno:
    def __init__(self, name, contig, strand, start, end, isoforms):
        self.name = name
        self.contig = contig
        self.strand = strand
        self.start = start
        self.end = end
        self.isoforms = isoforms
        self.len = abs(self.start-self.end)+1
        self.num_isoforms = len(self.isoforms)
        self.comp_mat = []
        self.region_len = []
        self.region_start = []
        self.region_end = []
    def addRegion(self, start, end, presence_isoform):



class RefSeqParser:
    def __init__(self, fname):
