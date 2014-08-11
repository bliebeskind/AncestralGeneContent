
def skip(src,line,dups=True):
	'''If dups is True, skip to duplications info, else skip to losses.'''
	if dups:
		key = "Duplication"
	else:
		key = "Species"
	while line:
		if not line.startswith(key):
			line = src.readline()
			continue
		else:
			line = src.readline()
			line = src.readline()
			break
	return src,line

def parse_dups(src,line,D):
	'''Add gains to dup/loss events dictionary.'''
	while line and line.strip() != '':
		line = line.strip().split(' ')
		line_list = [i for i in line if i != '']
		assert len(line_list) == 3
		node = line_list[1]
		if node not in D:
			D[node] = [1,0]
		else:
			D[node][0] += 1
		line = src.readline()
	return src,line,D
	
def parse_losses(src,line,D):
	'''Add losses to dup/loss events dictionary. Assumes each node is only
	represented once in losses section of Notung file.'''
	while line and line.strip() != '':
		line = line.strip().split(' ')
		line_list = [i for i in line if i != '']
		assert len(line_list) == 2
		node = line_list[0]
		num_losses = -int(line_list[1])
		if node not in D:
			D[node] = [0,num_losses]
		else:
			D[node][1] = num_losses
		line = src.readline()
	return src,line,D
	
def parse(infile):
	'''
	Parse Notung info file. Creates a dictionary with nodes as keys and a
	two member list as values. List is [<# of duplications>, <# of
	losses>]. Number of losses is a negative number.
	'''
	events = {}
	with open(infile) as f:
		line = f.readline()
		f,line = skip(f,line)
		f,line,events = parse_dups(f,line,events)
		f,line = skip(f,line,dups=False)
		f,line,events = parse_losses(f,line,events)
	return events