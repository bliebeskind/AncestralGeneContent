#! /usr/bin/env python

import sys, dendropy
import pandas as pd
import notung_info_parser

## info parser returnd a dictionary with gain/loss information for each node
## {'Chimp': [0, -1], # zero gains, 1 loss
## 'Homo': [0, 0],
## 'Mouse': [0, 0],
## 'Platypus': [0, 0],
## 'n2': [1, 0],
## 'n6': [0, 0],
## 'n8': [1, 0]}
##
## Currently needs to be modified to interact with class AncGenome better.
## Things to do are possibly remove tree construction functionality of 
## copy_number and batch_copy_number in favor of using functionality of class.
## Could also try making copy number into a generator, and generally look for
## chances to do this elsewhere.

def gl_branches(tree_file,info_file,format='newick'):
	'''Given a Notung gain/loss info file, annotate branches with 
	gains and losses in format: [gains, losses]. Losses will be 
	negative.'''
	tree = dendropy.Tree.get_from_path(tree_file,format)
	gl_D = notung_info_parser.parse(info_file)
	for node in tree.postorder_node_iter():
		if node.label:
			try:
				gl = gl_D[node.label]
			except KeyError:
				gl = "$" # placeholder for testing
			node.label = node.label + " %s" % str(gl)
		elif node.taxon:
			try:
				node.label = str(gl_D[str(node.taxon)])
			except KeyError:
				node.label = "$" # placeholder for testing
	return tree,gl_D
	
def get_name(node):
	'''Get name of node - either node label or taxon name.'''
	if node.label:
		return node.label
	elif node.taxon:
		return str(node.taxon)
	else:
		raise Exception("Node has neither label nor taxon")
		
def has_descendents(node,gl_D):
	'''If node has descendents that are in gl_D, return True, else
	return False'''
	for d in node.leaf_iter():
		if str(d.taxon) in gl_D:
			return True
	else:
		return False
	
def copy_number(tree_file,info_file,format='newick'):
	'''
	Annotate nodes with their gene copy number and return a tuple where
	the first item is the annotated (species) tree and the second is 
	mapping of nodes to counts
	
	Requires a species tree and the info output file from Notung. Uses
	information on duplications and losses to count up from the first 
	node that occured in the reconciled gene tree.
	'''
	tree = dendropy.Tree.get_from_path(tree_file,format)
	if not tree.is_rooted:
		tree.is_rooted == True
	gl_D = notung_info_parser.parse(info_file)
	node_counts = {}
	is_first = True # first node found that's in gl_D
	for node in tree.preorder_node_iter():
		name = get_name(node)
		if name in gl_D:
			if is_first:
				count = 1 + gl_D[name][0] + gl_D[name][1]
				node_counts[name] = count
				node.label = str(count)
				is_first = False
			else:
				count = int(node.parent_node.label) \
				+ gl_D[name][0] + gl_D[name][1]
				node_counts[name] = count
				node.label = str(count)
		else: # node not in info file, hence not in gene tree
			if node.parent_node == None or node.is_leaf(): # node is root or leaf
				node_counts[name] = 0
				node.label = "0"
			else:
				if has_descendents(node,gl_D):
					count = node.parent_node.label
					node_counts[node.label] = str(count)
					node.label = count
				else:
					node_counts[name] = 0
					node.label = "0"
	return tree, node_counts
	
def batch_copy_num(info_files,species_tree,write_trees=False,format='newick'):
	count_D = {}
	for f in info_files:
		gene = f.strip().split("_")[0]
		sys.stderr.write("Processing %s\n" % gene)
		try:
			annot_tree,node_counts = copy_number(species_tree,f,format)
			count_D[gene] = node_counts
		except IOError:
			raise Exception("Species tree or info_file not found")
		if write_trees:
			outfile = gene + "_Stree.nhx"
			annot_tree.write_to_path(outfile,'newick',
				suppress_leaf_node_labels=False)
	return pd.DataFrame(count_D)
			
