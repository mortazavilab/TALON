import pandas as pd
import argparse
import numpy as np
import csv

def get_args():

    desc = ('Replaces the starts or ends of transcripts in a GTF with the'
            ' longest alternatives (similar to the GENCODE model of'
            'calling transcripts)')
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-gtf', dest='gtf',
        help='TALON GTF to serve as the template to modify')
    parser.add_argument('-read_annot', dest='annot',
        help='Read annot file from TALON to extract raw read ends from')
    parser.add_argument("--datasets", "--d",  dest = "datasets_file",
        help = """A file indicating which datasets should be
                  included (one dataset name per line). Default is to include
                  all datasets.""", default='all')
    parser.add_argument('--mode', dest='mode', default='tss',
        help="Modify TSSs or TESs, 'tss', 'tes', 'both'. Default: 'tss'",
        choices={'tss', 'tes', 'both'})
    parser.add_argument('--novelty', dest='novelty',
        help="Whether to only modify ends from novel or all models, "+\
        "'all' or 'novel'. Default: 'all'", choices={'all', 'novel'}, default='all')
    parser.add_argument('-outprefix', '-o', dest='outprefix',
        help='Prefix for output file', default='talon')
    parser.add_argument('--verbose', '-v', action='store_true',
        default=False, help="Display in progress output")

    args = parser.parse_args()

    return args

# annot: TALON read annotation file
# how: 'tes' or 'tss' for calling ends or starts respectively
# novelty: 'all' or 'novel' based on which transcript models
#          you want to modify the ends of
# datasets: string file path of file with datasets to use when
#           calling longest ends from reads
def get_longest_ends(annot, how='tes', novelty='novel', datasets='a'):
    df = pd.read_csv(annot, sep='\t')

    if datasets != 'all':
        dataset_df = pd.read_csv(datasets, header=None, names=['dataset'])
        dataset_list = dataset_df['dataset'].tolist()
        for d in dataset_list:
            if d not in df.dataset.unique().tolist():
                raise ValueError("Dataset name {} not found in read_annot".format(d))

    if novelty == 'novel':
        df = df.loc[df.transcript_novelty != 'Known']

    fwd = df.loc[df.strand == '+']
    rev = df.loc[df.strand == '-']

    # furthest downstream for tes
    # if + strand, max coord of read end
    # if - strand, min coord of read end
    if how == 'tes':
        fwd = fwd[['transcript_ID', 'read_end']]
        fwd = fwd.groupby('transcript_ID').max().reset_index()
        rev = rev[['transcript_ID', 'read_end']]
        rev = rev.groupby('transcript_ID').min().reset_index()


    # furthest upstream for tss:
    # if + strand, min coord of read start
    # if - strand, max coord of read start
    elif how == 'tss':
        fwd = fwd[['transcript_ID', 'read_start']]
        fwd = fwd.groupby('transcript_ID').min().reset_index()
        rev = rev[['transcript_ID', 'read_start']]
        rev = rev.groupby('transcript_ID').max().reset_index()

    # concat fwd and rev
    df = pd.concat([fwd, rev])

    return df

# get the longest ends from the read annotation file
# annot: TALON read annotation file path
# how: 'tss' or 'tes', tss will find start ends and tes will find stop ends
# gtf: gtf file location
# ends: df with transcript_ID, end coordinate
# how: 'tss' or 'tes'
# opref: output file prefix
# verbose: display processing progress
# test: print out dataframe before and after editing
def replace_gtf_end_coords(gtf, ends, opref, how='tes', test=False, verbose=False):

    # read preexisting GTF and
    gtf_df = pd.read_csv(gtf, sep='\t', header=None, \
                names=['chr', 'source', 'entry_type', \
                       'start', 'stop', 'score', 'strand',\
                       'frame', 'fields'], comment='#')

    # get relevant values from fields
    gtf_df['transcript_id'] = np.nan
    gtf_df.loc[gtf_df.entry_type!='gene', 'transcript_id'] = gtf_df.loc[gtf_df.entry_type!='gene'].fields.str.split(pat='talon_transcript "', n=1, expand=True)[1]
    gtf_df.loc[gtf_df.entry_type!='gene', 'transcript_id'] = gtf_df.loc[gtf_df.entry_type!='gene'].transcript_id.str.split(pat='"', n=1, expand=True)[0]


    if how == 'tes':
        ends.columns = ['transcript_id', 'tes']
    elif how == 'tss':
        ends.columns = ['transcript_id', 'tss']

    # merge gtf_df with end information
#     ends.transcript_id = ends.transcript_id.astype('str')
    df = gtf_df.loc[gtf_df.transcript_id.notnull()]
    ends.transcript_id = ends.transcript_id.astype('str')
    gtf_df.transcript_id = gtf_df.transcript_id.astype('str')
    gtf_df = gtf_df.merge(ends, how='left', on='transcript_id')
    df.transcript_id = df.transcript_id.astype('str')
    df = df.merge(ends, how='inner')

    if test:
        print('Before editing')
        print(gtf_df[['transcript_id', 'entry_type', 'strand', 'start', 'stop', how]])

    # swap out read starts or ends for the longest ones
    tids = df.transcript_id.unique()
    for t, tid in enumerate(tids):
        if t % 1000 == 0 and verbose:
            print('Processing transcript {} of {}'.format(t, len(tids)))

        # fwd: swap out transcript "stop" and last exon "stop"
        # rev: swap out transcript "start" and last exon "start"
        if how == 'tes':
            # tes fwd
            ind = gtf_df.loc[(gtf_df.strand=='+')&(gtf_df.transcript_id==tid)].index.tolist()
            if ind:
                # stop of transcript for fwd
                i = ind[0]
                gtf_df.loc[i, 'stop'] = gtf_df.loc[i, 'tes']
                # stop of last exon for fwd
                i = ind[-1]
                gtf_df.loc[i, 'stop'] = gtf_df.loc[i, 'tes']

            # tes rev
            ind = gtf_df.loc[(gtf_df.strand=='-')&(gtf_df.transcript_id==tid)].index.tolist()
            if ind:
                # start of trancscript for rev
                i = ind[0]
                gtf_df.loc[i, 'start'] = gtf_df.loc[i, 'tes']
                # start of last exon for rev
                i = ind[-1]
                gtf_df.loc[i, 'start'] = gtf_df.loc[i, 'tes']
        # fwd: swap out transcript "start" and first exon "start"
        # rev: swap out transcript "stop" and first exon "stop"
        elif how == 'tss':
            # tss fwd
            ind = gtf_df.loc[(gtf_df.strand=='+')&(gtf_df.transcript_id==tid)].index.tolist()
            if ind:
                # start of transcript for fwd
                i = ind[0]
                gtf_df.loc[i, 'start'] = gtf_df.loc[i, 'tss']
                # start of first exon for fwd
                i = ind[1]
                gtf_df.loc[i, 'start'] = gtf_df.loc[i, 'tss']
            # tss rev
            ind = gtf_df.loc[(gtf_df.strand=='-')&(gtf_df.transcript_id==tid)].index.tolist()
            if ind:
                # stop of transcript for rev
                i = ind[0]
                gtf_df.loc[i, 'stop'] = gtf_df.loc[i, 'tss']
                # stop of first exon for rev
                i = ind[1]
                gtf_df.loc[i, 'stop'] = gtf_df.loc[i, 'tss']

    # now fix gene coordinates
    gtf_df['gene_id'] = gtf_df.loc[gtf_df.entry_type!='gene'].fields.str.split(pat='talon_gene "', n=1, expand=True)[1]
    gtf_df['gene_id'] = gtf_df.loc[gtf_df.entry_type!='gene'].gene_id.str.split(pat='"', n=1, expand=True)[0]

    # tes
    if how == 'tes':
        # fwd: replace "stop" of the gene with the maximum of the "stops"
        gene_ind = gtf_df.loc[(gtf_df.strand == '+')&(gtf_df.entry_type=='gene')&(gtf_df.tes.notnull())].index.tolist()
        if gene_ind:
            fwd = gtf_df.loc[(gtf_df.strand == '+')&(gtf_df.entry_type=='transcript')]
            if test:
                print('fwd')
                print(gtf_df.loc[gene_ind])
            gtf_df.loc[gene_ind, 'stop'] = gtf_df.loc[gene_ind].apply(lambda x: \
                fwd.loc[fwd.gene_id==x.gene_id, 'stop'].max(), axis=1)

        # rev: replace "start" of the gene with the minimum of the "starts"
        gene_ind = gtf_df.loc[(gtf_df.strand == '-')&(gtf_df.entry_type=='gene')&(gtf_df.tes.notnull())].index.tolist()
        if gene_ind:
            rev = gtf_df.loc[(gtf_df.strand == '-')&(gtf_df.entry_type=='transcript')]
            if test:
                print('rev')
                print(gtf_df.loc[gene_ind])
            gtf_df.loc[gene_ind, 'start'] = gtf_df.loc[gene_ind].apply(lambda x: \
                rev.loc[rev.gene_id==x.gene_id, 'start'].min(), axis=1)

    # tss
    elif how == 'tss':
        # fwd: replace "start" of the gene with the minimum of the "starts"
        gene_ind = gtf_df.loc[(gtf_df.strand == '+')&(gtf_df.entry_type=='gene')&(gtf_df.tss.notnull())].index.tolist()
        if gene_ind:
            fwd = gtf_df.loc[(gtf_df.strand == '+')&(gtf_df.entry_type=='transcript')]
            gtf_df.loc[gene_ind, 'start'] = gtf_df.loc[gene_ind].apply(lambda x: \
                fwd.loc[fwd.gene_id==x.gene_id, 'start'].min(), axis=1)

        # rev: replace "stop" of the gene with the maximum of the "stops"
        gene_ind = gtf_df.loc[(gtf_df.strand == '-')&(gtf_df.entry_type=='gene')&(gtf_df.tss.notnull())].index.tolist()
        if gene_ind:
            rev = gtf_df.loc[(gtf_df.strand == '-')&(gtf_df.entry_type=='transcript')]
            gtf_df.loc[gene_ind, 'stop'] = gtf_df.loc[gene_ind].apply(lambda x: \
                rev.loc[rev.gene_ind==x.gene_ind, 'stop'].max(), axis=1)

    if test:
        print()
        print('After editing')
        print(gtf_df[['transcript_id', 'entry_type', 'strand', 'start', 'stop', how]])

    cols=['chr', 'source', 'entry_type', \
          'start', 'stop', 'score', 'strand',\
           'frame', 'fields']
    gtf_df = gtf_df[cols]
    gtf_df['start'] = gtf_df['start'].astype('int')
    gtf_df['stop'] = gtf_df['stop'].astype('int')
    if test:
        fname = '{}_revised_{}_test.gtf'.format(opref, how)
    else:
        fname = '{}_revised_{}.gtf'.format(opref, how)
    gtf_df.to_csv(fname, sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)
    return gtf_df, fname

def main():

    args = get_args()
    gtf = args.gtf
    annot = args.annot
    mode = args.mode
    novelty = args.novelty
    oprefix = args.outprefix
    verbose = args.verbose
    datasets = args.datasets_file

    # first, call ends from the read annot file
    if mode == 'both':

        # tss first
        ends = get_longest_ends(annot, how='tss', novelty=novelty, datasets=datasets)
        df, fname = replace_gtf_end_coords(gtf, ends, oprefix,
            how='tss', verbose=verbose)

        # tes
        ends = get_longest_ends(annot, how='tes', novelty=novelty, datasets=datasets)
        df, fname = replace_gtf_end_coords(fname, ends, oprefix,
            how='tes', verbose=verbose)

    else:
        ends = get_longest_ends(annot, how=mode, novelty=novelty, datasets=datasets)
        df, fname = replace_gtf_end_coords(gtf, ends, oprefix,
            how=mode, verbose=verbose)
