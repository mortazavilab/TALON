import argparse
import os
import sqlite3

import numpy as np
import pandas as pd


def get_args():
    desc = (
        "Extracts the locations, novelty, and transcript assignments of"
        " exons/introns in a TALON database or GTF file. All positions "
        "are 1-based."
    )
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("--gtf", dest="gtf", default=None, help="TALON GTF file from which to extract exons/introns")
    parser.add_argument("--db", dest="db", default=None, help="TALON database from which to extract exons/introns")
    parser.add_argument(
        "--ref", dest="ref_gtf", help=("GTF reference file (ie GENCODE). Will be used to " "label novelty.")
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        help=(
            "Choices are 'intron' or 'exon' (default is 'intron'). "
            "Determines whether to include introns or exons in the "
            "output"
        ),
        default="intron",
    )
    parser.add_argument("--outprefix", dest="outprefix", help="Prefix for output file")

    args = parser.parse_args()

    if args.gtf and args.db:
        raise Exception("only input gtf or db")

    return args


# creates a dictionary of the last field of a gtf
# adapted from Dana Wyman
def get_fields(tab_fields):
    attributes = {}

    # remove trailing newline and split by semicolon
    description = tab_fields[-1].strip("\n")
    description = description.split(";")

    # Parse description
    for fields in description:
        if fields == "" or fields == " ":
            continue
        fields = fields.split()
        if fields[0] == "":
            fields = fields[1:]

        key = fields[0].replace('"', "")
        val = " ".join(fields[1:]).replace('"', "")

        attributes[key] = val

    # Put in placeholders for important attributes (such as gene_id) if they
    # are absent
    if "gene_id" not in attributes:
        attributes["gene_id"] = "NULL"

    return attributes


# create loc_df (for nodes), edge_df (for edges), and t_df (for paths)
def create_dfs_db(db):
    # make sure file exists
    if not os.path.exists(db):
        raise Exception("TALON db file not found. Check path.")

    # open db connection
    conn = sqlite3.connect(db)
    c = conn.cursor()

    # loc_df
    q = "SELECT loc.* FROM location loc"

    c.execute(q)
    locs = c.fetchall()

    loc_df = pd.DataFrame(locs, columns=["location_ID", "genome_build", "chrom", "position"])

    # do some df reformatting, add strand
    loc_df.drop("genome_build", axis=1, inplace=True)
    loc_df.rename({"location_ID": "vertex_id", "position": "coord"}, inplace=True, axis=1)
    loc_df.vertex_id = loc_df.vertex_id.map(int)

    # edge_df
    q = """SELECT * FROM edge """

    c.execute(q)
    edges = c.fetchall()

    edge_df = pd.DataFrame(edges, columns=["edge_id", "v1", "v2", "edge_type", "strand"])
    edge_df.v1 = edge_df.v1.map(int)
    edge_df.v2 = edge_df.v2.map(int)
    edge_df["talon_edge_id"] = edge_df.edge_id
    edge_df["edge_id"] = edge_df.apply(lambda x: (int(x.v1), int(x.v2)), axis=1)

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

    t_df["tid"] = np.asarray(tids)
    t_df["path"] = np.asarray(paths)

    t_df = create_dupe_index(t_df, "tid")
    t_df = set_dupe_index(t_df, "tid")

    # furnish the last bit of info in each df
    t_df["path"] = [[int(n) for n in path] for path in get_db_vertex_paths(paths, edge_df)]
    loc_df = create_dupe_index(loc_df, "vertex_id")
    loc_df = set_dupe_index(loc_df, "vertex_id")

    edge_df.drop("talon_edge_id", axis=1, inplace=True)
    edge_df = create_dupe_index(edge_df, "edge_id")
    edge_df = set_dupe_index(edge_df, "edge_id")

    return loc_df, edge_df, t_df

    # create loc_df (nodes), edge_df (edges), and t_df (transcripts) from gtf
    # adapted from Dana Wyman and TALON


def create_dfs_gtf(gtf_file):
    # make sure file exists
    if not os.path.exists(gtf_file):
        raise Exception("GTF file not found. Check path.")

    # depending on the strand, determine the stard and stop
    # coords of an intron or exon
    def find_edge_start_stop(v1, v2, strand):
        if strand == "-":
            start = max([v1, v2])
            stop = min([v1, v2])
        elif strand == "+":
            start = min([v1, v2])
            stop = max([v1, v2])
        return start, stop

    # dictionaries to hold unique edges and transcripts
    transcripts = {}
    exons = {}

    with open(gtf_file) as gtf:
        for line in gtf:
            # ignore header lines
            if line.startswith("#"):
                continue

            # split each entry
            line = line.strip().split("\t")

            # get some fields from gtf that we care about
            chrom = line[0]
            entry_type = line[2]
            start = int(line[3])
            stop = int(line[4])
            strand = line[6]
            fields = line[-1]

            # transcript entry
            if entry_type == "transcript":
                attributes = get_fields(line)
                tid = attributes["transcript_id"]
                gid = attributes["gene_id"]

                # add transcript to dictionary
                transcript = {tid: {"gid": gid, "tid": tid, "strand": strand, "exons": []}}
                transcripts.update(transcript)

            # exon entry
            elif entry_type == "exon":
                attributes = get_fields(line)
                start, stop = find_edge_start_stop(start, stop, strand)
                eid = "{}_{}_{}_{}_exon".format(chrom, start, stop, strand)
                tid = attributes["transcript_id"]

                # add novel exon to dictionary
                if eid not in exons:
                    edge = {eid: {"eid": eid, "chrom": chrom, "v1": start, "v2": stop, "strand": strand}}
                    exons.update(edge)

                    # add this exon to the transcript's list of exons
                if tid in transcripts:
                    transcripts[tid]["exons"].append(eid)

    # once we have all transcripts, make loc_df
    locs = {}
    vertex_id = 0
    for edge_id, edge in exons.items():
        chrom = edge["chrom"]
        strand = edge["strand"]

        v1 = edge["v1"]
        v2 = edge["v2"]

        # exon start
        key = (chrom, v1)
        if key not in locs:
            locs[key] = vertex_id
            vertex_id += 1
        # exon end
        key = (chrom, v2)
        if key not in locs:
            locs[key] = vertex_id
            vertex_id += 1

    # add locs-indexed path to transcripts, and populate edges
    edges = {}
    for _, t in transcripts.items():
        t["path"] = []
        strand = t["strand"]
        t_exons = t["exons"]

        for i, exon_id in enumerate(t_exons):
            # pull some information from exon dict
            exon = exons[exon_id]
            chrom = exon["chrom"]
            v1 = exon["v1"]
            v2 = exon["v2"]
            strand = exon["strand"]

            # add current exon and subsequent intron
            # (if not the last exon) for each exon to edges
            key = (chrom, v1, v2, strand)
            v1_key = (chrom, v1)
            v2_key = (chrom, v2)
            edge_id = (locs[v1_key], locs[v2_key])
            if key not in edges:
                edges[key] = {"edge_id": edge_id, "edge_type": "exon"}

            # add exon locs to path
            t["path"] += list(edge_id)

            # if this isn't the last exon, we also needa add an intron
            # this consists of v2 of the prev exon and v1 of the next exon
            if i < len(t_exons) - 1:
                next_exon = exons[t_exons[i + 1]]
                v1 = next_exon["v1"]
                key = (chrom, v2, v1, strand)
                v1_key = (chrom, v1)
                edge_id = (locs[v2_key], locs[v1_key])
                if key not in edges:
                    edges[key] = {"edge_id": edge_id, "edge_type": "intron"}

    # turn transcripts, edges, and locs into dataframes
    locs = [{"chrom": key[0], "coord": key[1], "vertex_id": vertex_id} for key, vertex_id in locs.items()]
    loc_df = pd.DataFrame(locs)

    edges = [
        {
            "v1": item["edge_id"][0],
            "v2": item["edge_id"][1],
            "strand": key[3],
            "edge_id": item["edge_id"],
            "edge_type": item["edge_type"],
        }
        for key, item in edges.items()
    ]
    edge_df = pd.DataFrame(edges)

    transcripts = [{"tid": key, "gid": item["gid"], "path": item["path"]} for key, item in transcripts.items()]
    t_df = pd.DataFrame(transcripts)

    # final df formatting
    loc_df = create_dupe_index(loc_df, "vertex_id")
    loc_df = set_dupe_index(loc_df, "vertex_id")
    edge_df = create_dupe_index(edge_df, "edge_id")
    edge_df = set_dupe_index(edge_df, "edge_id")
    t_df = create_dupe_index(t_df, "tid")
    t_df = set_dupe_index(t_df, "tid")

    return loc_df, edge_df, t_df


# convert talon query into edge path
def get_db_edge_paths(paths):
    edge_paths = []
    for p in paths:
        if p[1] == None:
            edge_paths.append([p[0]])
        else:
            edge_paths.append([p[0], *[int(i) for i in p[1].split(",")], p[2]])
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
            else:
                path.append(entry.v2.values[0])
        vertex_paths.append(path)
    return vertex_paths


# creates the duplicate index
def create_dupe_index(df, ind_name):
    df[ind_name + "_back"] = df[ind_name]
    return df


def add_coord_info(edge_df, loc_df):
    edge_df["chrom"] = edge_df.apply(lambda x: loc_df.loc[x.v1, "chrom"], axis=1)
    edge_df["start"] = edge_df.apply(lambda x: loc_df.loc[x.v1, "coord"], axis=1)
    edge_df["stop"] = edge_df.apply(lambda x: loc_df.loc[x.v2, "coord"], axis=1)

    return edge_df


def subset_edges(edge_df, mode="intron"):
    sjs = edge_df[edge_df.apply(lambda x: True if x.edge_type == mode else False, axis=1)]
    return sjs


def determine_sj_novelty(ref_edge_df, edge_df):
    # Merge known starts from ref_edge_df with the query edges
    ref_edge_df["start_known"] = True
    edge_df = edge_df.merge(
        ref_edge_df[["chrom", "start", "strand", "start_known"]], how="left", on=["chrom", "strand", "start"]
    )
    edge_df.fillna(value=False, inplace=True)

    # Merge known ends from ref_edge_df with the query edges
    ref_edge_df["stop_known"] = True
    edge_df = edge_df.merge(
        ref_edge_df[["chrom", "stop", "strand", "stop_known"]], how="left", on=["chrom", "strand", "stop"]
    )
    edge_df.fillna(value=False, inplace=True)

    # Now determine whether the edge in whole has been seen before
    ref_edge_df["combination_known"] = True
    edge_df = edge_df.merge(
        ref_edge_df[["chrom", "start", "stop", "strand", "combination_known"]],
        how="left",
        on=["chrom", "strand", "start", "stop"],
    )
    edge_df.fillna(value=False, inplace=True)

    return edge_df


# renames old index dupe column in df and resets the index
def reset_dupe_index(df, ind_name):
    df.rename({ind_name: ind_name + "_back"}, inplace=True, axis=1)
    df.reset_index(inplace=True)
    return df


# set index, rename dupe index in df
def set_dupe_index(df, ind_name):
    df.set_index(ind_name, inplace=True)
    df.rename({ind_name + "_back": ind_name}, inplace=True, axis=1)
    return df


def format_edge_df(edge_df):
    edge_df.reset_index(drop=True, inplace=True)
    edge_df.drop(["edge_type", "v1", "v2"], axis=1, inplace=True)
    return edge_df


def find_tids_from_sj(edge_df, t_df, mode="intron"):
    if mode == "exon":
        t_df["edges"] = t_df.apply(lambda x: [(x.path[i], x.path[i + 1]) for i in range(len(x.path[:-1]))][::2], axis=1)
    elif mode == "intron":
        t_df["edges"] = t_df.apply(
            lambda x: [(x.path[i], x.path[i + 1]) for i in range(len(x.path[:-1]))][1::2], axis=1
        )
    edge_df["tids"] = edge_df.apply(lambda x: add_tids_to_sj(x, t_df), axis=1)
    edge_df.reset_index(drop=True, inplace=True)
    edge_df.drop("edge_id", inplace=True, axis=1)

    return edge_df


def add_tids_to_sj(x, t_df):
    return ",".join([tid for tid, edges in zip(t_df.tid, t_df.edges) if x.edge_id in edges])


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
    edge_df = determine_sj_novelty(ref_edge_df, edge_df)
    edge_df = find_tids_from_sj(edge_df, t_df, mode=args.mode)

    edge_df = edge_df.rename(columns={"tids": "transcript_ids"})
    edge_df.to_csv(
        "{}_{}s.tsv".format(args.outprefix, args.mode),
        sep="\t",
        index=False,
        columns=[
            "chrom",
            "start",
            "stop",
            "strand",
            "start_known",
            "stop_known",
            "combination_known",
            "transcript_ids",
        ],
    )


if __name__ == "__main__":
    main()
