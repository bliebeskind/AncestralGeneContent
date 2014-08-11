#! /usr/bin/env

import pandas as pd
import dendropy
import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties

def node_path(taxon,tree):
	taxon_node = tree.find_node_with_taxon_label(taxon)
	path = []
	for node in taxon_node.ancestor_iter():
		if node.label and node.parent_node != None:
			path.append(node.label)
		elif node.taxon:
			path.append(str(node.taxon))
		elif node.parent_node == None:
			continue
		else:
			raise Exception("Tree has unlabeled nodes")
	path.reverse()
	return path

def path_df(path,df):
	path_df = pd.DataFrame(columns=df.columns)
	for i in path:
		path_df = path_df.append(df.ix[i])
	return path_df
	
def print_paths(tree_file,gene_counts,format='newick'):
	tree = dendropy.Tree.get_from_path(tree_file,format)
	tree.is_rooted = True
	taxa = [str(i.taxon) for i in tree.leaf_iter()]
	all_counts = pd.read_csv(gene_counts,index_col=0)
	for i in taxa:
		print "Processing %s" % i
		path = node_path(i,tree)
		df = path_df(path,all_counts)
		plot = df.plot(legend=False)
		plot.legend(df.columns,ncol=3,loc='best')
		plt.savefig(i+"_path.pdf")