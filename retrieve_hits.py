#! /usr/bin/env python

import os,sys
import cPickle as pickle
from Bio import SeqIO

## This needs to be fixed, it's way too ugly now. Find a better (not ad hoc)
## way of supplying the path to the fasta files. Also maybe make it work as a
## stream output instead of building a SeqRecord list.

def get_pickled_hits(infile):
	'''Load in dictionary made with hits_pickler'''
	with open(infile) as f:
		return pickle.load(f)

def walk(start="."):
	'''Return list of tuples from os.walk'''
	return [i for i in os.walk(start)]
		
def path_to_file(name,walk):
	'''Returns list where sole entry is path to file that starts with "name"'''
	pathlist = []
	for dirpath, dirnames, filenames in walk:
		for f in filenames:
			if f.startswith(name) and not f.endswith(".txt"):
				pathlist.append(os.path.join(dirpath,f))
				assert len(pathlist) < 2,\
				"Found more than one file matching name: %s" % name
	assert len(pathlist) > 0, "Didn't find any files matching '%s'" % name
	sys.stderr.write("Retrieving from %s\n" % pathlist[0])
	return pathlist
		
def get_hits(path, hits):
	'''Return list of fasta SeqIO objects. These are specified from dictionary 
	object "hits," which has dictionary entries where keys are sequence ids. 
	Gets these from file in "path"'''
	hits = hits.keys() # get list of hits
	outrecs = []
	records = SeqIO.parse(path,'fasta')
	for rec in records:
		if rec.id in hits:
			outrecs.append(rec)
	assert len(outrecs) == len(hits), "Didn't find all the hits"
	return outrecs
		
def hits_from_files(D,path_to_fastas=None):
	'''Make and return list of fasta SeqIO objects from dictionary D made with 
	hits_pickler.'''
	all_hits = []
	for key in D.keys():
		if path_to_fastas:
			pathwalk = walk(path_to_fastas)
			filepath = path_to_file(key,pathwalk)
			new_hits = get_hits(filepath[0],D[key])
		else:
			filepath = path_to_file(key,pathwalk)
			new_hits = get_hits(filepath[0],D[key])
		all_hits += new_hits
	return sorted(all_hits)
	
if __name__ == '__main__':
	infile = sys.argv[1]
	outfile = sys.argv[2]
	try:
		path_to_fastas = sys.argv[3]
	except IndexError:
		path_to_fastas = None
	pickled_hits = get_pickled_hits(infile)
	species = 0
	hits = 0
	for k in pickled_hits.keys():
		species +=1
		hits += len(pickled_hits[k].keys())
	sys.stderr.write("%i total species\n" % species)
	sys.stderr.write("%i total hits\n" % hits)
	if path_to_fastas:
		all_hits = hits_from_files(pickled_hits,path_to_fastas)
	assert len(all_hits) == hits
	sys.stderr.write("Found all the hits!\n")
	SeqIO.write(all_hits,outfile,'fasta')
