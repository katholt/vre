##########################
#dr.mark.schultz@gmail.com
#github: schultzm
#130215
#Happy Friday the 13th!
##########################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import os
import argparse

#example execute "python parse.py -g AAAK01000287.gbk AE016832.gbk > name_of_log_file.txt"
#set up the arguments parser to deal with the command line input
parser = argparse.ArgumentParser(description = "Python/Biopython script for outputting the DNA sequences of all 'genes' in a Genbank input to separate files. Input format Genbank only.  Output format, any biopython supported format (e.g., Genbank, fasta, phylip, nexus).")
parser.add_argument('-g', '--genbank', nargs='+',required=True, help = 'Input names of Genbank file(s), separated by whitespace, to pull the productname and sequence slice from.')
parser.add_argument('-o', '--out_format', default='fasta', required=False, help = "Output file format. Default is 'fasta'.")
args = parser.parse_args()

#parse the genbank file and write every 'gene' in the file out to a separate file with filename and header = Accesion number + gene name + gene location
def extract_genes(file, outformat):
	with open(file, 'r') as handle:
		print "\n********"
		print "Parsing",file,"..."
		for record in SeqIO.parse(handle, "genbank", generic_dna):
			for i in range(0, len(record.features)):
				if record.features[i].type == "gene":
					#the values stored in the feature variable are needed to get the seq slices
					feature = record.features[i]
					gene_name = record.features[i].qualifiers["gene"]
					#get rid of spaces in name
					gene_name = gene_name[0].replace(" ","_")
					accession_number = str(os.path.splitext(file)[0])
					new_record_id = accession_number+"_"+gene_name
					my_seq_location = record.features[i].location
					#for examples on extracting seqs from genbank see http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/#features
					#thanks Peter for the blog
					my_seq_string = str(feature.extract(record.seq))
					record_new = SeqRecord(Seq(my_seq_string, generic_dna), id=new_record_id, description="| 'Original seq coordinates in refseq "+accession_number+" were: "+str(my_seq_location)+"'")
					outfile_name = new_record_id+"."+outformat
					print "Extracting",gene_name,"gene in location",str(my_seq_location),"from GenBank record",new_record_id,"and writing to", outfile_name
					with open(outfile_name, 'w') as output_handle:
						SeqIO.write(record_new, outfile_name, outformat)

print "\nAttempting to process files: "+str(args.genbank)

#execute the extract_genes function on every file in the input list
for i in args.genbank:
	extract_genes(i, args.out_format)
