from talon import talon_label_reads as tlr
import pandas as pd

def test_pool_outputs():
    """ Given some SAM files and some log files, check the concatenation
        process. """
    indir = "talon_label_reads/test_inputs/pool_test"
    outprefix = "scratch/pool_test"

    tlr.pool_outputs(indir, outprefix)

    sam = outprefix + "_labeled.sam"
    log = outprefix + "_read_labels.tsv"

    # Check content: SAM file should have 2 header lines and 6 reads
    n_header = 0
    sam_entries = 0
    read_ids = set()
    expected_read_ids = set(["read_1", "read_2", "read_3", "read_4", "read_5",
                             "read_6"])
    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                n_header += 1
            else:
                sam_entries += 1
                read_ids.add(line.split()[0])
    assert n_header == 2
    assert sam_entries == 6
    assert read_ids == expected_read_ids
    
    # Check content: Log should have a header line and 6 entries
    log_data = pd.read_csv(log, sep="\t", header = 0)
    assert len(log_data) == 6
    assert set(list(log_data.read_name)) == expected_read_ids
