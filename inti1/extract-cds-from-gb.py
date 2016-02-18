import sys

#Script to take a file of proteins in GenBank/GenPept format, examine
#their annotation, and use this to download their CDS from the NCBI.
#Written and tested on Biopython 1.49, on 2008/01/19
# by Peter biopython at maubp.freeserve.co.uk
# Modified by Adina Howe, 2016

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
#Edit this next line (and read the NCBI Entrez usage guidelines)
Entrez.email = "adina@iastate.edu"

def get_nuc_by_name(name, start=None, end=None) :
   """Fetches the sequence from the NCBI given an identifier name.

   Note start and end should be given using one based counting!

   Returns a Seq object."""
   record = SeqIO.read(Entrez.efetch("nucleotide",
                                     id=name.strip(),
                                     seq_start=start,
                                     seq_stop=end,
                                     retmode="text",
                                     rettype="fasta"), "fasta")
   return record.seq

def get_nuc_from_coded_by_string(source) :
   """Fetches the sequence from the NCBI for a "coded_by" string.

   e.g. "NM_010510.1:21..569" or "AF376133.1:<1..>553"
   or "join(AB061020.1:1..184,AB061020.1:300..1300)"
   or "complement(NC_001713.1:67323..68795)"

   Note - joins and complements are handled by recursion.

   Returns a Seq object."""

   if source.startswith("complement(") :
       assert source.endswith(")")
       #For simplicity this works by recursion
       return get_nuc_from_coded_by_string(source[11:-1]).reverse_complement()

   if source.startswith("join(") :
       assert source.endswith(")")
       #For simplicity this works by recursion.
       #Note that the Seq object (currently) does not have a join
       #method, so convert to strings and join them, then go back
       #to a Seq object:
       return Seq("".join(str(get_nuc_from_coded_by_string(s)) \
                          for s in source[5:-1].split(",")))

   if "(" in source or ")" in source \
   or source.count(":") != 1 or source.count("..") != 1 :
       raise ValueError("Don't understand %s" % repr(source))
   name, loc = source.split(":")
   #Remove and ignore any leading < or > for fuzzy locations which
   #indicate the full CDS extends beyond the region sequenced.
   start, end = [int(x.lstrip("<>")) for x in loc.split("..")]
   #We could now download the full sequence, and crop it locally:
   #return get_nuc_by_name(name)[start-1:end]
   #However, we can ask the NCBI to crop it and then download
   #just the bit we need!
   return get_nuc_by_name(name,start,end)

def get_nuc_record(protein_record, table="Standard") :
   if not isinstance(protein_record, SeqRecord) :
       raise TypeError("Expect a SeqRecord as the protein_record.")
   feature = None
   #print dir(protein_record)
   desc = protein_record.description
   name = protein_record.id
   prot =  protein_record.seq
   #print "----"
   for f in protein_record.features :
       if f.type == "CDS" and "coded_by" in f.qualifiers :
           feature = f
           break
   if feature :
       #This is the good situation, there is a precise "coded_by" string
       #Check this CDS feature is for the whole protein:
       assert feature.location.start.position == 0
       assert feature.location.end.position == len(protein_record)
       source = feature.qualifiers["coded_by"][0]
       #print "Using %s" % source
       #print feature
       #return SeqRecord(Seq(''))
       nuc = get_nuc_from_coded_by_string(source)
       #See if this included the stop codon - they don't always!
       if str(nuc[-3:].translate(table)) == "*" :
           nuc = nuc[:-3]
       return name, source, nuc, prot, desc




file_i = SeqIO.parse(open(sys.argv[1]), "genbank")
for p in file_i:
   id, id2, nuc_seq, prot_seq, desc = get_nuc_record(p, table="Standard")
   file_op = open(str(id) + '.pep.fa', 'w')
   file_on = open(str(id) + '.nuc.fa', 'w')
   file_on.write(">%s\t%s\n" % (id, desc))
   file_on.write("%s\n" % nuc_seq)
   file_op.write(">%s\t%s\n" % (id, id2))
   file_op.write("%s\n" % prot_seq)
