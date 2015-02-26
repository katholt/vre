import string, re
import os, sys, subprocess
import collections
from optparse import OptionParser
from Bio import Phylo

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file (newick)", default="")
	parser.add_option("-n", "--names", action="store", dest="namesfile", help="names file (csv)", default="")
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file (newick)", default="")
	parser.add_option("-a", "--add", action="store", dest="add", help="create names by adding the two columns (rather than replacing with column 2 text)", default="0")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
			
	if options.tree=="":
		DoError("No tree file specified (-t)")	
		
	if options.namesfile=="":
		DoError("No names file specified (comma-sep, column1 = ID, column2 = new names) (-n)")	
		
	if options.outputfile=="":
		DoError("No output file specified (newick format tree with new names) (-n)")	
		
	idtable = {}
	f = file(options.namesfile, "r")
	for line in f:
		fields = line.rstrip().split(",")
		if options.add=="1":
			idtable[fields[0]] = fields[0] + "_" + fields[1]
		else:
			idtable[fields[0]] = fields[1]
		
	t = Phylo.read(options.tree, "newick") # read in newick tree
	for node in t.get_terminals():
		name = node.name
  		if name in idtable:
  			node.name = idtable[name]
  		else:
  			node.name = name
  		print "replaced node " + name + " with " + node.name
  	Phylo.write(t,options.outputfile,"newick")