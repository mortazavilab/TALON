# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Fairlie Reese
# -----------------------------------------------------------------------------
# create_anndata_from_database.py is a utility that outputs the abundance
# for each transcript in the TALON database across datasets in AnnData format.

import sqlite3
import itertools
import operator
from optparse import OptionParser
from pathlib import Path
import scanpy
import numpy as np
import pandas as pd
import anndata
from scipy.sparse import csr_matrix

from . import filter_talon_transcripts as filt
from .. import dstruct as dstruct
from .. import length_utils as lu
from . import post_utils as putils
from . import ab_utils as autils
from .. import query_utils as qutils
from .. import talon as talon


def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")

    parser.add_option("--annot", "-a", dest = "annot",
        help = """Which annotation version to use. Will determine which
                  annotation transcripts are considered known or novel
                  relative to. Note: must be in the TALON database.""",
        type = "string")

    parser.add_option("--pass_list", dest = "pass_list",
                      help = "Pass list file of transcripts to include in the \
                              output. First column should be TALON gene ID, \
                              second column should be TALON transcript ID",
                      metavar = "FILE", type = "string", default = None)

    parser.add_option("--build", "-b", dest = "build",
        help = "Genome build to use. Note: must be in the TALON database.",
        type = "string")

    parser.add_option('--gene', dest='gene_level',
        help='Output AnnData on the gene level rather than the transcript',
        action='store_true')

    parser.add_option("--datasets", "-d",  dest = "dataset_file",
        help = """Optional: A file indicating which datasets should be
                  included (one dataset name per line). Default is to include
                  all datasets.""",
        metavar = "FILE", type = "string", default = None)

    parser.add_option("--o", dest = "ofile", help = "Output file name",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def assign_novelties(df, d, order, how):
    """
    Assign novelty types based on a priority order
    for either genes or transcripts

    Parameters:
        df (pandas DataFrame): Long-form DataFrame of
            gene / transcript statuses for each entry
        d (dict): Dictionary mapping novelty category names
            to attribute / value combinations in df
        order (list of str): List of novelty category names
            to assign to a gene / transcript in priority order
        how (str): {'gene', 'transcript'}

    Returns:
        df (pandas DataFrame): DataFrame indexed by gene /
            transcript ID with novelty information
    """
    if how == 'gene':
        nov_col = 'gene_novelty'
        cols = [nov_col]
    elif how == 'transcript':
        nov_col = 'transcript_novelty'
        cols = [nov_col, 'ISM_subtype']

    # assign gene or transcript novelty
    df = df.pivot(index='ID', columns=['attribute'], values=['value'])
    df = df.droplevel(0, axis=1)
    df.columns.name = ''

    for key, value in d.items():
        df[key] = False

        # in cases where we're filtering out a lot,
        # not all novelty types will be represented
        if value[0] in df.columns:
            df.loc[df[value[0]]==value[1], key] = True
            df.drop(value[0], axis=1, inplace=True)

    df[nov_col] = np.nan
    for o in order:
        df.loc[(df[nov_col].isnull())&(df[o]==True), nov_col] = o

    # assign ism subtype if needed
    if how == 'transcript':
        df['ISM_subtype'] = np.nan
        df.loc[(df.ISM_subtype.isnull())&(df['ISM-prefix'])&(df['ISM-suffix']), 'ISM_subtype'] = 'Both'
        df.loc[(df.ISM_subtype.isnull())&(df['ISM-prefix']), 'ISM_subtype'] = 'Prefix'
        df.loc[(df.ISM_subtype.isnull())&(df['ISM-suffix']), 'ISM_subtype'] = 'Suffix'
        df.loc[df.ISM_subtype.isnull(), 'ISM_subtype'] = 'None'

    # reduce cols
    df = df[cols]
    df.reset_index(inplace=True)

    return df

def get_transcript_novs(db, tids):
    """
    Get transcript novelties and ISM subtypes from a TALON db

    Parameters:
        db (str): Path to TALON db
        tids (list of int): List of transcript IDs to include

    Returns:
        df (pandas DataFrame): DF with novelties and
            ISM subtypes from a TALON db
    """

    # attributes to search for
    nov_col_dict = {'Known': ('transcript_status', 'KNOWN'),
                    'ISM': ('ISM_transcript', 'TRUE'),
                    'ISM-prefix': ('ISM-prefix_transcript', 'TRUE'),
                    'ISM-suffix': ('ISM-suffix_transcript', 'TRUE'),
                    'NIC': ('NIC_transcript', 'TRUE'),
                    'NNC': ('NNC_transcript', 'TRUE'),
                    'Antisense': ('antisense_transcript', 'TRUE'),
                    'Intergenic': ('intergenic_transcript', 'TRUE'),
                    'Genomic': ('genomic_transcript', 'TRUE')}
    order = ['ISM', 'NIC', 'NNC', 'Antisense',
             'Intergenic', 'Genomic', 'Known']
    attr_list = [val[0] for key, val in nov_col_dict.items()]
    attrs = qutils.format_for_IN(attr_list)

    # transcripts to search for
    transcript_query = qutils.format_for_IN(tids)

    with sqlite3.connect(db) as conn:
        query = f"""SELECT ID, attribute, value
                    FROM transcript_annotations
                    WHERE attribute IN {attrs}
                    AND ID IN {transcript_query}
                 """
        df = pd.read_sql_query(query, conn)

    df = assign_novelties(df, nov_col_dict, order, 'transcript')

    return df

def get_gene_novs(db, gids):
    """
    Get gene novelties from a TALON db

    Parameters:
        db (str): Path to TALON db
        gids (list of int): Gene IDs to include

    Returns:
        df (pandas DataFrame): DF with novelties from a TALON db
    """

    # attributes to search for
    nov_col_dict = {'Known': ('gene_status', 'KNOWN'),
                   'Intergenic': ('intergenic_novel', 'TRUE'),
                   'Antisense': ('antisense_gene', 'TRUE')}
    order = ['Antisense', 'Intergenic', 'Known']
    attr_list = [val[0] for key, val in nov_col_dict.items()]
    attrs = qutils.format_for_IN(attr_list)

    # genes to search for
    gene_query = qutils.format_for_IN(gids)

    with sqlite3.connect(db) as conn:
        query = f"""SELECT ID, attribute, value
                    FROM gene_annotations
                    WHERE attribute IN {attrs}
                    AND ID IN {gene_query}
                 """
        df = pd.read_sql_query(query, conn)
    df = assign_novelties(df, nov_col_dict, order, 'gene')

    return df

def get_g_t_names(db, annot, tids):
    """
    Get names / IDs of genes / transcripts from TALON db

    Parameters:
        db (str): Path to TALON db
        annot (str): Name of annotation in TALON db
        tids (list of int): List of transcript IDs to include

    Returns:
        df (pandas DataFrame): DataFrame holding name / ID
            info for each gene / transcript
    """
    t_query = qutils.format_for_IN(tids)

    # get information that we want for each transcript, stuff that
    # would be output in the abundance table
    with sqlite3.connect(db) as conn:
        query = f"""
            SELECT
                t.gene_ID,
                t.transcript_ID,
                ga_id.value AS annot_gene_id,
                ta_id.value AS annot_transcript_id,
                ga_name.value AS annot_gene_name,
                ta_name.value AS annot_transcript_name,
                t.n_exons
            FROM transcripts t
                LEFT JOIN gene_annotations ga_id ON t.gene_ID = ga_id.ID
                    AND ga_id.annot_name = '{annot}'
                    AND ga_id.attribute = 'gene_id'
                LEFT JOIN transcript_annotations ta_id ON t.transcript_ID = ta_id.ID
                    AND ta_id.annot_name = '{annot}'
                    AND ta_id.attribute = 'transcript_id'
                LEFT JOIN gene_annotations ga_name ON t.gene_ID = ga_name.ID
                    AND ga_name.annot_name = '{annot}'
                        AND ga_name.attribute = 'gene_name'
                LEFT JOIN transcript_annotations ta_name ON t.transcript_ID = ta_name.ID
                    AND ta_name.annot_name = '{annot}'
                        AND ta_name.attribute = 'transcript_name'
                WHERE t.transcript_ID in {t_query}
            """
        df = pd.read_sql_query(query, conn)

    return df

def get_var_info(db, annot, build, tids=None, gids=None, gene_level=False):
    """
    Get info about names, IDs, novelty categories, etc. for each gene
    and transcript in a talon DB

    Parameters:
        db (str): Path to TALON db
        annot (str): Name of annotation in TALON db
        build (str): Name of genome build in TALON db
        tids (list of int): Internal transcript IDs to include
        gids (list of int): Internal gene IDs to include
        gene_level (bool): Whether to return info on the gene level

    Returns:
        df (pandas DataFrame): DataFrame with metadata about
            each gene and transcript in a TALON db
    """

    # get names / ids of transcripts / genes
    df = get_g_t_names(db, annot, tids)
    prefix = autils.fetch_naming_prefix(db)
    n_places = autils.fetch_n_places(db)

    # make names for novel genes / transcripts
    # determine how many missing digits there are
    # repeat '0' for that many spaces for each gene / transcript ID
    df['zero'] = '0'
    df['n_gid_zero_to_add'] = n_places-df.gene_ID.astype(str).str.len()
    df['temp_gid'] = prefix+'G'+df['zero'].str.repeat(df['n_gid_zero_to_add'])+df['gene_ID'].astype(str)
    df['n_tid_zero_to_add'] = n_places-df.transcript_ID.astype(str).str.len()
    df['temp_tid'] = prefix+'T'+df['zero'].str.repeat(df['n_tid_zero_to_add'])+df['transcript_ID'].astype(str)

    df['temp'] = df.temp_gid.str.len()
    if len(df['temp'].unique().tolist()) != 1:
        raise ValueError('Problem naming genes')
    df['temp'] = df.temp_tid.str.len()
    if len(df['temp'].unique().tolist()) != 1:
        raise ValueError('Problem naming transcripts')

    # drop extra stuff
    drop_cols = ['zero', 'n_gid_zero_to_add',
                 'n_tid_zero_to_add', 'temp']
    df.drop(drop_cols, axis=1, inplace=True)

    # # add gene / transcript names / ids
    # df[['temp_gid', 'temp_tid']] = df.apply(lambda x: talon.construct_names(x.gene_ID,
    #                      x.transcript_ID,
    #                      prefix,
    #                      n_places),
    #                  axis=1, result_type='expand')

    # replace null gene names / ids
    inds = df.loc[df.annot_gene_id.isnull()].index
    df.loc[inds, 'annot_gene_id'] = df.loc[inds, 'temp_gid']
    inds = df.loc[df.annot_gene_name.isnull()].index
    df.loc[inds, 'annot_gene_name'] = df.loc[inds, 'temp_gid']

    # replace null transcript names / ids
    inds = df.loc[df.annot_transcript_id.isnull()].index
    df.loc[inds, 'annot_transcript_id'] = df.loc[inds, 'temp_tid']
    inds = df.loc[df.annot_transcript_name.isnull()].index
    df.loc[inds, 'annot_transcript_name'] = df.loc[inds, 'temp_tid']

    # remove temp cols
    df.drop(['temp_gid', 'temp_tid'], axis=1, inplace=True)

    # add transcript len
    t_lens = pd.DataFrame.from_dict(autils.get_transcript_lengths(db, build),
                                                orient='index',
                                                columns=['length'])
    df = df.merge(t_lens, how='left', left_on='transcript_ID', right_index=True)

    # add gene novelty
    g_df = get_gene_novs(db, gids)
    df = df.merge(g_df, how='left', left_on='gene_ID', right_on='ID')
    df.drop('ID', axis=1, inplace=True)

    # add transcript novelty / ism subtype
    t_df = get_transcript_novs(db, tids)
    df = df.merge(t_df, how='left', left_on='transcript_ID', right_on='ID')
    df.drop('ID', axis=1, inplace=True)

    # column order
    order = ['gene_ID', 'transcript_ID', 'annot_gene_id',
             'annot_transcript_id', 'annot_gene_name',
             'annot_transcript_name', 'n_exons', 'length',
             'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    df = df[order]

    # gene level -- drop columns that are only relevant to transcripts
    # and drop duplicated entries
    if gene_level:
        drop_cols = ['transcript_ID', 'annot_transcript_id',
                     'annot_transcript_name', 'length',
                     'transcript_novelty', 'ISM_subtype',
                     'n_exons']
        df.drop(drop_cols, axis=1, inplace=True)
        df.drop_duplicates(inplace=True)

    return df

def get_obs_info(db, dataset_file):
    """
    Get metadata table for each dataset in TALON

    Parameters:
        db (str): Path to TALON db
        dataset_file (str): Path to file with datasets
            to include

    Returns:
        df (pandas DataFrame): DataFrame with info
            about each dataset queried for
    """
    datasets = autils.fetch_dataset_list(dataset_file, db)
    datasets_query = qutils.format_for_IN(datasets)
    with sqlite3.connect(db) as conn:
        query = f"""
                 SELECT dataset_name, sample, platform
                 FROM dataset WHERE dataset_name IN {datasets_query}
                 """
        df = pd.read_sql_query(query, conn)
        df.rename({'dataset_name': 'dataset'}, axis=1, inplace=True)
    return df

def get_X_info(db, obs, var, gene_level=False):
    """
    Get sparse matrix representation of gene or transcript counts
    from the TALON db

    Parameters:
        db (str): Path to TALON db
        obs (pandas DataFrame): Pandas DataFrame with information about each
            dataset / sample to include
        var (pandas DataFrame): Pandas DataFrame with information about each
            gene or transcript to include
        gene_level (bool): Whether to compute counts on the gene / transcript level
    """

    dataset_str = qutils.format_for_IN(obs.dataset.unique().tolist())

    # filter on genes
    if gene_level:
        var_col = 'gene_ID'
        feat_str = qutils.format_for_IN(var[var_col].unique().tolist())
        query = f"""SELECT t.gene_ID, ab.transcript_ID, ab.dataset, ab.count
                    FROM abundance as ab
                    LEFT JOIN transcripts as t
                        ON t.transcript_ID = ab.transcript_ID
                    WHERE t.gene_ID in {feat_str}
                    AND ab.dataset in {dataset_str}
                 """

    # filter on transcripts
    else:
        var_col = 'transcript_ID'
        feat_str = qutils.format_for_IN(var[var_col].unique().tolist())
        query = f"""SELECT transcript_ID, dataset, count
                    FROM abundance WHERE transcript_ID in {feat_str}
                    AND dataset in {dataset_str}
                 """
    with sqlite3.connect(db) as conn:
        df = pd.read_sql_query(query, conn)

    # # remove genes / transcripts w/ 0 counts
    # df = df.loc[df['count'] > 0]

    # sum over transcripts from the same gene / dataset
    if gene_level:
        df.drop('transcript_ID', axis=1, inplace=True)
        df = df.groupby(['gene_ID', 'dataset']).sum().reset_index()

    # make categories based on ordering of obs and var tables
    obs_col = 'dataset'
    obs_cat = pd.api.types.CategoricalDtype(obs[obs_col], ordered=True)
    if obs_cat.categories.tolist() != obs[obs_col].tolist():
        raise ValueError('Problem with dataset names')
    var_cat = pd.api.types.CategoricalDtype(var[var_col], ordered=True)
    if var_cat.categories.tolist() != var[var_col].tolist():
        raise ValueError('Problem with feature IDs')

    # create sparse matrix representation without
    # inflating
    row = df[obs_col].astype(obs_cat).cat.codes
    col = df[var_col].astype(var_cat).cat.codes
    X = csr_matrix((df['count'], (row, col)), \
                   shape=(obs_cat.categories.size,
                          var_cat.categories.size))

    # # code to inflate matrix
    # dfs = pd.SparseDataFrame(X, \
    #                      index=obs_cat.categories, \
    #                      columns=var_cat.categories, \
    #                      default_fill_value=0)

    return X

def main():
    options = getOptions()
    db = options.database
    annot = options.annot
    build = options.build

    pass_list_file = options.pass_list
    dataset_file = options.dataset_file
    ofile = options.ofile
    gene_level = options.gene_level

    # make sure that the input database exists
    if not Path(db).exists():
        raise ValueError("Database file '%s' does not exist!" % db)

    autils.check_annot_validity(annot, db)
    autils.check_build_validity(build, db)

    # determine which transcripts to include
    pass_list = putils.handle_filtering(db,
                                        annot,
                                        True,
                                        pass_list_file,
                                        dataset_file)
    gids = [i[0] for i in list(set(pass_list))]
    tids = [i[1] for i in list(set(pass_list))]

    # get obs, var, and X tables
    var = get_var_info(db, annot, build, tids, gids, gene_level)
    obs = get_obs_info(db, dataset_file)
    X = get_X_info(db, obs, var, gene_level)

    # assemble adata
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.write(ofile)

if __name__ == '__main__':
    main()
