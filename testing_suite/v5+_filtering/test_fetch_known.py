import make_minimal_db_for_filtering as mmdb

def test_init():
    """ Initialize the minimal database needed to run the next tests.
        For testing inclusion of Knowns
         1. A known transcript (id = 1) (dataset 1) with low fracA
         2. A known transcript (id = 2) (dataset 2) with high fracA
        For testing selective exclusion of genomic transcripts and min_dataset
        and min_count criteria
         3. A genomic transcript (id = 3) (dataset 1) with low fracA
         4. A genomic transcript (id = 3) (dataset 2) with low fracA
         5. A genomic transcript (id = 3) (dataset 3) with low fracA

    """

    db_file = "test.db" # TODO: change path

    # Add reads. Fields that are not relevant for this purpose are set to None
    known = [(1, 1, 1, "read_1", "dataset_1", None, None, None, None, None, None,
              None, 0.2, None, None),
             (2, 1, 2, "read_2", "dataset_2", None, None, None, None, None, None,
              None, 0.7, None, None)]
    genomic = [(3, 1, 3, "read_3", "dataset_1", None, None, None, None, None, None,
                None, 0.2, None, None),
               (4, 1, 3, "read_4", "dataset_2", None, None, None, None, None, None,
               None, 0.2, None, None),
               (5, 1, 3, "read_5", "dataset_3", None, None, None, None, None, None,
              None, 0.2, None, None)]

    reads = known + genomic

    # Datasets
    datasets = [(1, "dataset_1", "test", "test"),
                (2, "dataset_2", "test", "test"),
                (3, "dataset_3", "test", "test")] 

    # Annotations
    annotations = [(1, "", "", "transcript_status", "KNOWN"),
                   (2, "", "", "transcript_status", "KNOWN"),
                   (3, "", "", "transcript_status", "NOVEL"),
                   (3, "", "", "genomic_transcript", "TRUE")]

    mmdb.make_minimal_db_for_filtering(db_file, reads, datasets, annotations)


