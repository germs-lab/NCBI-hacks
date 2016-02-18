import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


genome=SeqIO.read(sys.argv[1], 'genbank')


l = []
n = 0
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    org = record.annotations["source"]
    for feat in genome.features:
        if feat.type == "CDS":
            print feat.coded_by
            start = feat.location.start.position
            end = feat.location.end.position
            pos = [start, end]
            print feat
            l.append(pos)
            print '>' + sys.argv[1].split('.')[0] + org + ' '+ '16S rRNA gene' + str(n)                

            print feat.extract(genome.seq)
            n = n + 1

