import os
import pandas as pd
import glob
import seaborn as sns
import argparse
import time
import sqlite3
import numpy as np


def get_args():

	desc = ('Extracts the locations, novelty, and transcript assignments of'
                ' exons/introns in a TALON database or GTF file. All positions '
                'are 1-based.')
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('--gtf', dest='gtf', default=None,
		help = 'TALON GTF file from which to extract exons/introns')
	parser.add_argument('--db', dest='db', default=None,
		help = 'TALON database from which to extract exons/introns')
	parser.add_argument('--ref', dest='ref_gtf', 
		help = ('GTF reference file (ie GENCODE). Will be used to '
                        'label novelty.'))
	parser.add_argument('--mode', dest='mode', 
		help= ("Choices are 'intron' or 'exon' (default is 'intron'). "
			"Determines whether to include introns or exons in the "
			"output"), default='intron')
	parser.add_argument('--outprefix', dest='outprefix',
		help = 'Prefix for output file')


	args = parser.parse_args()

	if args.gtf and args.db: 
		raise Exception('only input gtf or db')

	return args


# get value associated with keyword in the 9th column of gtf
def get_field_value(key, fields):
	if key not in fields:
		return None
	else:
		return fields.split(key+' "')[1].split()[0].replace('";','')

# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
def create_dfs_db(db):

	# make sure file exists
	if not os.path.exists(db):
		raise Exception('TALON db file not found. Check path.')

	# open db connection
	conn = sqlite3.connect(db)
	c = conn.cursor()

	# loc_df
	q = 'SELECT loc.* FROM location loc'

	c.execute(q)
	locs = c.fetchall()

	loc_df = pd.DataFrame(locs,
		columns=['location_ID', 'genome_build',
				 'chrom', 'position'])

	# do some df reformatting, add strand
	loc_df.drop('genome_build', axis=1, inplace=True)
	loc_df.rename({'location_ID': 'vertex_id',
				   'position': 'coord'},
				   inplace=True, axis=1)
	loc_df.vertex_id = loc_df.vertex_id.map(int)

	# edge_df
	q = """SELECT * FROM edge """

	c.execute(q)
	edges = c.fetchall()

	edge_df = pd.DataFrame(edges, 
		columns=['edge_id', 'v1', 'v2',
				 'edge_type', 'strand'])
	edge_df.v1 = edge_df.v1.map(int)
	edge_df.v2 = edge_df.v2.map(int)
	edge_df['talon_edge_id'] = edge_df.edge_id
	edge_df['edge_id'] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

	# t_df
	t_df = pd.DataFrame()

	# get tid, gid, gname, and paths
	q = """SELECT ga.value, ta.value,
				  t.start_exon, t.jn_path, t.end_exon,
				  t.start_vertex, t.end_vertex
			FROM gene_annotations ga 
			JOIN transcripts t ON ga.ID=t.gene_ID
			JOIN transcript_annotations ta ON t.transcript_ID=ta.ID
			WHERE ta.attribute='transcript_id'
			AND (ga.attribute='gene_name' 
			OR ga.attribute='gene_id')
		"""

	c.execute(q)
	data = c.fetchall()

	# get fields from each transcript and add to dataframe
	gids, tids, paths = zip(*[(i[0], i[1], i[2:]) for i in data[::2]])
	gnames = [i[0] for i in data[1::2]]
	paths = get_db_edge_paths(paths)

	t_df['tid'] = np.asarray(tids)
	t_df['gid'] = np.asarray(gids)
	t_df['gname'] = np.asarray(gnames)
	t_df['path'] = np.asarray(paths)

	t_df = create_dupe_index(t_df, 'tid')
	t_df = set_dupe_index(t_df, 'tid')

	# furnish the last bit of info in each df
	t_df['path'] = [[int(n) for n in path]
					 for path in get_db_vertex_paths(paths, edge_df)]
	loc_df['strand'] = loc_df.apply(lambda x:
			 get_db_strand(x, edge_df), axis=1)
	loc_df = create_dupe_index(loc_df, 'vertex_id')
	loc_df = set_dupe_index(loc_df, 'vertex_id')

	edge_df.drop('talon_edge_id', axis=1, inplace=True)
	edge_df = create_dupe_index(edge_df, 'edge_id')
	edge_df = set_dupe_index(edge_df, 'edge_id')

	return loc_df, edge_df, t_df

# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
def create_dfs_gtf(gtf):

	# make sure file exists
	if not os.path.exists(gtf):
		raise Exception('GTF file not found. Check path.')

	# get dfs by parsing through gtf
	loc_df = pd.DataFrame(columns=['chrom', 'coord', 'strand','vertex_id'])
	loc_df.set_index(['chrom', 'coord'], inplace=True)

	edge_df = pd.DataFrame(columns=['edge_id', 'edge_type',
								    'strand', 'v1', 'v2'])
	t_df = pd.DataFrame(columns=['tid', 'gid',
								 'gname', 'path'])

	# loop initialization
	vertex_id = 0
	transcript_paths = []
	transcript_path = []

	with open(gtf, 'r') as infile:
		for line in infile:

			# skip header lines
			if '##' in line: continue

			line = line.strip().split('\t')

			# gene entry 
			if line[2] == 'gene':
				curr_gid = get_field_value('gene_id', line[-1])
				curr_gname = get_field_value('gene_name', line[-1])

			# transcript entry
			elif line[2] == 'transcript':
				curr_tid = get_field_value('transcript_id', line[-1])
				
				# start a new transcript path
				if transcript_path != []:

					# add to list of transcript paths and transcript df 
					transcript_paths.append(transcript_path)
					t_df = t_df.append({'tid': prev_tid,
								 'gid': prev_gid,
								 'gname': prev_gname,
								 'path': transcript_path},
								 ignore_index=True)

				transcript_path = []

				# reset some stuff
				terminal_loc = True
				exon = 0
				intron = 1

			# exon entry
			elif line[2] == 'exon':

				# get exon info 
				curr_chr = line[0]
				curr_start = line[3]
				curr_stop = line[4]
				curr_strand = line[6]
				
				if curr_strand == '+': coords = [curr_start, curr_stop]
				else: coords = [curr_stop, curr_start]
				
				for c in coords:

					ind = (curr_chr, int(c))

					# loc not in loc_df already
					if ind not in loc_df.index.tolist():

						attr = {'vertex_id': vertex_id,	   
								'coord': int(c),
								'strand': curr_strand, 'chrom': curr_chr}

						# update loc_df and increment vertex_id
						loc_df.reset_index(inplace=True)
						loc_df = loc_df.append(attr, ignore_index=True)
						loc_df.set_index(['chrom', 'coord'], inplace=True)

						curr_loc = int(vertex_id)
						vertex_id += 1

					# loc was already added to graph
					else: curr_loc = int(loc_df.loc[ind].vertex_id)	
	
					# add an edge to previous loc if not terminal loc 
					# and if the edge doesn't already exist
					if not terminal_loc:
						curr_edge = (prev_loc, curr_loc)
						
						if curr_edge not in edge_df.edge_id.to_list():
							attrs = {'edge_id': (curr_edge[0], curr_edge[1]),
								     'v1': curr_edge[0],
									 'v2': curr_edge[1], 
									 'strand': curr_strand}
							if exon: attrs.update({'edge_type': 'exon'})
							elif intron: attrs.update({'edge_type': 'intron'})

							edge_df = edge_df.append(attrs, ignore_index=True)

					# update transcript path with each loc 
					transcript_path.append(curr_loc)
					prev_loc = curr_loc
					prev_tid = curr_tid
					prev_gid = curr_gid
					prev_gname = curr_gname
					terminal_loc = False
					
					# exon or intron
					exon = abs(exon-1)
					intron = abs(intron-1)
					
	# append last transcript info
	transcript_paths.append(transcript_path)
	t_df = t_df.append({'tid': curr_tid,
					    'gid': curr_gid,
					    'gname': curr_gname,
						'path': transcript_path},
						ignore_index=True)

	# label node/edge types and finish formatting dfs correctly
	loc_df.reset_index(inplace=True)
	loc_df = create_dupe_index(loc_df, 'vertex_id')
	loc_df = set_dupe_index(loc_df, 'vertex_id')

	loc_df.coord = loc_df.coord.map(int)

	t_df = create_dupe_index(t_df, 'tid')
	t_df = set_dupe_index(t_df, 'tid')

	edge_df = create_dupe_index(edge_df, 'edge_id')
	edge_df = set_dupe_index(edge_df, 'edge_id')

	return loc_df, edge_df, t_df

# convert talon query into edge path
def get_db_edge_paths(paths):
	edge_paths = []
	for p in paths:
		if p[1] == None:
			edge_paths.append([p[0]])
		else:
			edge_paths.append(
				[p[0], *[int(i) for i in p[1].split(',')], p[2]])
	return edge_paths

# convert edge path to vertex path
def get_db_vertex_paths(paths, edge_df):
	vertex_paths = []
	for p in paths: 
		path = []
		for i, e in enumerate(p): 
			entry = edge_df.loc[edge_df.talon_edge_id == e]
			if i == 0:
				path.extend([entry.v1.values[0], entry.v2.values[0]])
			else: path.append(entry.v2.values[0])
		vertex_paths.append(path)
	return vertex_paths

# get the strand of each vertex
def get_db_strand(x, edge_df):
	# use v1 or v2 depending on where vertex is in edge
	try: 
		strand = edge_df.loc[edge_df.v1 == x.vertex_id, 'strand'].values[0]
	except:
		strand = edge_df.loc[edge_df.v2 == x.vertex_id, 'strand'].values[0]
	return strand

# creates the duplicate index
def create_dupe_index(df, ind_name):
	df[ind_name+'_back'] = df[ind_name]
	return df

def add_coord_info(edge_df, loc_df):
	edge_df['chrom'] = edge_df.apply(lambda x: loc_df.loc[x.v1, 'chrom'], axis=1)
	edge_df['start'] = edge_df.apply(lambda x: loc_df.loc[x.v1, 'coord'], axis=1)
	edge_df['stop'] = edge_df.apply(lambda x: loc_df.loc[x.v2, 'coord'], axis=1)

	return edge_df

def subset_edges(edge_df, mode='intron'):
	sjs = edge_df[edge_df.apply(
		lambda x: True if x.edge_type == mode else False, axis=1)]
	return sjs

def determine_sj_novelty(ref_loc_df, ref_edge_df, loc_df, edge_df):

	drop_cols = ['vertex_id']
	ref_loc_df.drop(drop_cols, axis=1, inplace=True)
	loc_df.drop(drop_cols, axis=1, inplace=True)

	ref_loc_df.reset_index(drop= True,inplace=True)
	ref_loc_df['annotated'] = True
	loc_df.reset_index(drop=True, inplace=True)

	merged_locs = ref_loc_df.merge(loc_df, how='right', 
		on=['chrom', 'strand', 'coord'])
	merged_locs.fillna(value=False, inplace=True)

	merged_locs.set_index(['chrom', 'strand', 'coord'], inplace=True)

	edge_df['start_known'] = edge_df.apply(
		lambda x: merged_locs.loc[(x.chrom, x.strand, x.start), 'annotated'],
		axis=1)
	edge_df['stop_known'] = edge_df.apply(
		lambda x: merged_locs.loc[(x.chrom, x.strand, x.stop), 'annotated'],
		axis=1)

        # Combined novelty (i.e. has start been seen with stop before)        
	ref_edge_df.reset_index(drop= True,inplace=True)
	ref_edge_df.drop('edge_id', axis = 1, inplace=True)
	ref_edge_df['combination_known'] = True
	merged_edge = ref_edge_df.merge(edge_df, how='right',
                on=['chrom', 'strand', 'start', 'stop'])
	merged_edge.fillna(value=False, inplace=True)

	#edge_df = create_dupe_index(edge_df, 'edge_id')
	#edge_df = set_dupe_index(edge_df, 'edge_id')

	return merged_edge

# renames old index dupe column in df and resets the index
def reset_dupe_index(df, ind_name):
	df.rename({ind_name: ind_name+'_back'}, inplace=True, axis=1)
	df.reset_index(inplace=True)
	return(df)

# set index, rename dupe index in df
def set_dupe_index(df, ind_name):
	df.set_index(ind_name, inplace=True)
	df.rename({ind_name+'_back': ind_name}, inplace=True, axis=1)
	return(df)

def format_edge_df(edge_df):
	edge_df.reset_index(drop=True, inplace=True)
	edge_df.drop(['edge_type', 'v1', 'v2'], axis=1, inplace=True)
	return edge_df

def find_tids_from_sj(edge_df, t_df, mode='intron'):
	if mode == 'exon':
		t_df['edges'] = t_df.apply(
			lambda x: [(x.path[i], x.path[i+1]) for i in range(len(x.path[:-1]))][::2],
			axis=1)
	elif mode == 'intron':
		t_df['edges'] = t_df.apply(
			lambda x: [(x.path[i], x.path[i+1]) for i in range(len(x.path[:-1]))][1::2],
			axis=1)
	edge_df['tids'] = edge_df.apply(lambda x: add_tids_to_sj(x, t_df), axis=1)
	edge_df.reset_index(drop=True, inplace=True)
	edge_df.drop('edge_id', inplace=True, axis=1)

	return edge_df

def add_tids_to_sj(x, t_df):
	# return [tid for tid, edges in zip(t_df.tid, t_df.edges) if x.edge_id in edges]
	return ','.join([tid for tid, edges in zip(t_df.tid, t_df.edges) if x.edge_id in edges])


def main():
	args = get_args()

	ref_loc_df, ref_edge_df, ref_t_df = create_dfs_gtf(args.ref_gtf)
	ref_edge_df = add_coord_info(ref_edge_df, ref_loc_df)
	ref_edge_df = subset_edges(ref_edge_df, mode=args.mode)
	ref_edge_df = format_edge_df(ref_edge_df)

	if args.db: 
		loc_df, edge_df, t_df = create_dfs_db(args.db)

	elif args.gtf:
		loc_df, edge_df, t_df = create_dfs_gtf(args.gtf)

	edge_df = add_coord_info(edge_df, loc_df)
	edge_df = subset_edges(edge_df, mode=args.mode)
	edge_df = format_edge_df(edge_df)
	edge_df = determine_sj_novelty(ref_loc_df, ref_edge_df, loc_df, edge_df)
	edge_df = find_tids_from_sj(edge_df, t_df, mode=args.mode)

	edge_df.to_csv('{}_{}_sj_summary.tsv'.format(args.outprefix, args.mode), 
			sep='\t', index=False)


if __name__ == '__main__':
	main()
