import evolmap
import notung_gl_tree
import pandas as pd
import dendropy
import copy

class AncGenome:
	'''
	Class for analyzing ancestral genome content data from EvolMap and Notung.
	'''
	
	def __init__(self):
		self.node_count = None 	# pandas DataFrame of node counts
		self.node_dist = None	# pandas DataFrame of node distribution
		self.tree = None		# dendropy Tree object
		
	def load_evolmap(self,file_list,tree_file,format='newick'):
		'''
		Load node counts from several EvolMap analyses. Populates 
		variable self.node_counts and self.tree
		
		file_list: A python list of ancestor_results files to load from
		tree_file: Tree file with labeled nodes. Must have all the taxa
			that are in the ancestor results file. Can specify 'format'
		'''
		self.node_count,self.tree = \
		evolmap.count_df(file_list,tree_file,False,format)
		
	def load_notung(self,info_files,tree_file,format='newick'):
		'''
		Load node counts from Notung .info files. Populates self.node_counts
		and self.tree.
		
		info_files: python list of .info files. Must start with: "GeneName_"
		tree_file: Tree file with labeled nodes. Can specify format.
		'''
		self.node_count, self.tree = \
		notung_gl_tree.batch_copy_num(info_files,tree_file,format)
		
	def normalize(self):
		'''
		Populate self.node_dist with a normalized DataFrame. Each row
		(current or extent genome) will sum to one and will therefore be 
		a distribution over gene types.
		'''
		T = self.node_count.T
		norm = (T/T.sum()).T
		self.node_dist = norm.fillna(0)
		
	def _calc_gl(self,count,parent_count,kind):
		'''Calculate gain, loss, or fold loss between parent node and node.'''
		if kind.lower() == 'gain':
			gain = count - parent_count
			if gain < 0:
				gain = 0.0
			return gain
		elif kind.lower() == 'loss':
			loss = parent_count - count
			if loss < 0:
				loss = 0.0
			return loss
		elif kind.lower() == 'foldloss':
			fold_loss = (parent_count - count) / parent_count
			if fold_loss > 1:
				raise Exception("Fold loss cannot be greater than 1")
			elif fold_loss < 0:
				fold_loss = 0
			return fold_loss
		else:
			raise Exception("kind %s not recognized" % kind)
		
		
	def gl_branches(self,gene,kind='gain'):
		'''
		Create and return dendropy Tree object whose branches are 
		proportional to gene gain or loss. Can choose three types:
			kind = "gain" or "loss" or "foldloss"
		'''
		tree = copy.deepcopy(self.tree)
		for node in tree.preorder_node_iter():
			if node.parent_node==None:	# EvolMap doesn't have counts for 
				node.edge_length = 0.0		# root, this may be a problem
			elif node.is_internal():
				assert node.label, "Interior node has no label"
				try:
					count = self.node_count[gene][node.label]
				except KeyError:
					raise Exception("Node %s not found in node counts" % (node.label))
				parent_count = float(self.node_count[gene][node.parent_node.label])
				node.edge_length = self._calc_gl(count,parent_count,kind)
			elif node.is_leaf():
				assert node.taxon, "Leaf has no associated taxon"
				try:
					count = self.node_count[gene][str(node.taxon)]
					parent_count = \
						float(self.node_count[gene][node.parent_node.label])
					node.edge_length = self._calc_gl(count,parent_count,kind)
				except KeyError:
					raise Exception("Node %s not found in node counts" 
										% (node.label))
		return tree
		
	def path_df(self,taxon,tree_file=None,format='newick'):
		'''Return a DataFrame of channel counts for the path from root
		to specified taxon. Can be plotted as a line graph to show the mode
		of duplication and loss.'''
		if tree_file:
			tree = dendropy.Tree.get_from_path(tree_file, format)
		else:
			tree = self.tree
		taxon_node = tree.find_node_with_taxon_label(taxon)
		path_df = pd.DataFrame(columns=self.node_count.columns)
		path = [i for i in taxon_node.ancestor_iter()]
		path.reverse()
		for node in path:
			if node.label and node.parent_node != None:
				path_df = path_df.append(self.node_count.ix[node.label])
			elif node.taxon:
				path_df = path_df.append(self.node_count.ix[str(node.taxon)])
			elif node.parent_node == None:
				continue
			else:
				raise Exception("Tree has unlabeled nodes")
		return path_df
		
	def rename_nodecount_index(self,rename_dict):
		'''Rename the indices of self.node_count using an input dictionary
		mapping old names to new names.'''
		new = self.node_count.T
		new = new.rename(columns = rename_dict)
		self.node_count = new.T
