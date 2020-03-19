import argparse
import pandas as pd

def get_args():

	desc = 'Fixes a GTF with no genes'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-gtf', '-g', dest='gtf',
		help='gtf to fix')
	args = parser.parse_args()

	return args

# check what entries are missing in the gtf
def is_bad_gtf(gtffile):

	missing_gene = False
	missing_trans = False

	# how many lines are useless lines
	with open(gtffile, 'r') as infile:
		for i, line in enumerate(infile): 
			if '##' not in line:
				break
	skiprows = [j for j in range(0, i)]

	df = pd.read_csv(gtffile, sep='\t', usecols=[2], skiprows=skiprows)
	categories = df.iloc[:,0].unique()

	# print(categories)

	# what are we missing?
	if 'gene' not in categories:
		missing_gene = True
	if 'transcript' not in categories: 
		missing_trans = True

	return (missing_gene, missing_trans)

# get value associated with keyword in the 9th column of gtf
def get_field_value(key, fields):
    if key not in fields:
        return None
    else:
        return fields.split(key+' "')[1].split()[0].replace('";','')

def construct_new_entry(prev_line, coords, entry_type):

	# print('Constructing new {} entry'.format(entry_type))
	
	# add gene or transcript type, coords, and len
	prev_line[2] = entry_type
	prev_line[3] = min(coords)
	prev_line[4] = max(coords)
	prev_line[7] = '.'

	# change the fields to reflect what type we are now
	new_fields = ''
	fields = prev_line[-1]
	gid = get_field_value('gene_id', fields)
	new_fields += 'gene_id "{}";'.format(gid)

	# if there's a gene name add it too
	gname = get_field_value('gene_name', fields)
	if gname:
		new_fields += 'gene_name "{}";'.format(gname)

	if entry_type == 'transcript':
		tid = get_field_value('transcript_id', fields)
		new_fields += ' transcript_id "{}";'.format(tid)
		
	prev_line[-1] = new_fields
	prev_line = format_to_write(prev_line)

	return prev_line

def make_ofile_name(matfile, prefix=None):
	fname = matfile.split('.gtf')[0]
	if prefix:
		fname += '_'
		fname += prefix
	fname += '_reformatted.gtf'
	return fname

def format_to_write(line): 
	return ''.join('\t'.join([str(i) for i in line])+'\n')

def main():

	args = get_args()
	gtffile = args.gtf

	(missing_gene, missing_transcript) = is_bad_gtf(gtffile)

	print('Missing transcript :  {}'.format(missing_transcript))

	# if nothing is missing, you good!
	if not missing_gene and not missing_transcript: 
		print('GTF has both gene and transcript entries. Nothing to add.')
		return

	# loop through this thing
	infile = open(gtffile, 'r')
	outfile = open(make_ofile_name(gtffile), 'w')

	curr_gid = ''
	curr_gid_coords = []

	curr_tid = ''
	curr_tid_coords = []

	first_transcript = True
	first_exon = True

	gene_list = []
	transcript_list = []

	prev_line = ''

	# relevant entries
	entries = ['exon']

	if missing_gene: 
		entries.append('transcript')

	if missing_gene or missing_transcript:

		for line in infile: 

			# skip the dumb header lines
			if line.startswith('#'):
				continue

			line = line.strip().split('\t')
			fields = line[-1]

			gid = get_field_value('gene_id', fields)
			tid = get_field_value('transcript_id', fields)

			if line[2] in entries:

				# set variables if first entry
				if first_exon: 
					curr_gid = gid
					curr_tid = tid

					curr_gid_coords = [int(line[3]), int(line[4])]
					curr_tid_coords = [int(line[3]), int(line[4])]

					first_exon = False

					prev_line = line

				# found a new transcript
				elif missing_transcript and tid != curr_tid: 

					# create transcript entry and dump to current gene list
					new_entry = construct_new_entry(
						prev_line, curr_tid_coords, 'transcript')
					transcript_list = new_entry+''.join(transcript_list) 

					gene_list += transcript_list 
					transcript_list = ''
					curr_tid_coords = []

				if missing_gene and gid != curr_gid:
					
					# create gene entry and write current gene list
					new_entry = construct_new_entry(
						prev_line, curr_gid_coords, 'gene')
					gene_list = new_entry+''.join(gene_list)

					gene_list += ''.join(transcript_list)
					transcript_list = ''
					curr_tid_coords = []

					outfile.write(gene_list)
					gene_list = ''
					curr_gid_coords = []

				# update loop vars
				curr_gid = gid
				curr_tid = tid
				curr_gid_coords.append(int(line[3]))
				curr_gid_coords.append(int(line[4]))
				curr_tid_coords.append(int(line[3]))
				curr_tid_coords.append(int(line[4]))

				prev_line = line

			# regardless, append to list of entries to write
			transcript_list += format_to_write(line)

		# if we've reached the end of the file
		# create transcript entry and dump to current gene list
		if missing_transcript:
			new_entry = construct_new_entry(
				prev_line, curr_tid_coords, 'transcript')
			transcript_list = new_entry+''.join(transcript_list)

		gene_list += transcript_list 
		transcript_list = ''

		# create gene entry and write current gene list
		if missing_gene:
			new_entry = construct_new_entry(
				prev_line, curr_gid_coords, 'gene')

		gene_list = new_entry+''.join(gene_list)
		outfile.write(gene_list)
		gene_list = ''


		infile.close()
		outfile.close()

if __name__ == '__main__':
	main()