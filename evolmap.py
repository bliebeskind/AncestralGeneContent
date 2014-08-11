#! /usr/bin/env

import dendropy
import pandas as pd

def parse_anc_results(anc_results):
	'''Return generator of lines in ancestor_results file. Yields dictionaries
	holding values with fields as keys.'''
	with open(anc_results) as f:
		is_first = True
		for line in f:
			if is_first:
				is_first = False
				line = line.strip().split("\t")
				fields = line
				continue
			line = line.strip().split("\t")
			line_D = {}
			for pos,value in enumerate(line):
				line_D[fields[pos]] = value
			if line_D["Ancestor name"] == '':
				break
			yield line_D
			
def get_counts(anc_results,tree_file,mapping=None,format='newick'):
	'''
	Parse ancestor_results file from Evolmap and return generator of
	nodes mapped to counts. User can supply a dictionary mapping the strings
	in the first columns of the ancestor_results file to node names. If
	mapping is not supplied, a tree file with node names must be supplied.
	This tree file must have all the taxa that are in ancestor_results file.
	'''
	tree = dendropy.Tree.get_from_path(tree_file,format)
	yield tree
	for line in parse_anc_results(anc_results):   
		taxa = line["Ancestor name"].strip().split("_")
		if len(taxa) > 1: # is internal
			tree_node = tree.mrca(taxon_labels=taxa)
			assert tree_node.label, \
				"Node label not found. Taxa: %s" % str(taxa)
			node = tree_node.label
			try:
				yield node, line["Present loci"]
			except KeyError: # no loci found in anc_results file
				if tree_node.parent_node == None: 	# OK if root
					yield node, line['Sym-bets']	# yield sym-bets
				else: # But all other nodes should have counts
					raise Exception("No loci found for taxa: %s\n" % str(taxa))
		else: # is leaf
			try:
				yield taxa[0], line["Present loci"]
			except KeyError:
				raise Exception("No loci found for taxon: %s\n" % taxa[0])
                   
def count_df(file_list,tree_file,normalized=False,format='newick'):
	'''Given a list of EvolMap output files of format Gene1.ancestor_results,
	and a tree with node labels, return a pandas DataFrame of gene counts for
	each node and leaf and the dendropy Tree object'''
	is_first = True
	D = {}
	for f in file_list:
		gene = f.split(".")[0]
		count_gen = get_counts(f,tree_file,mapping=None)
		tree = count_gen.next() # first value is tree
		counts = pd.Series({node:count for node,count in count_gen})
		D[gene] = counts
	df = pd.DataFrame(D).astype("float")
	if normalized:
		T = df.T
		dist_df = (T/T.sum()).T
		return dist_df.fillna(0),tree
	else:
		return df,tree
