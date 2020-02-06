import pyfaidx
from talon import talon_label_reads as tlr

def test_firstbase_plus_strand():
    """ Fetch the first base of this genome on the + strand.
            ACTGACTGACTGAAATAAGAAACTGACTG 
        The correct answer is A
    """
    genome_file = "talon_label_reads/test_inputs/toy_genome.fa"
    genome = pyfaidx.Fasta(genome_file, sequence_always_upper=True,
                           one_based_attributes=False)
    assert tlr.fetch_seq("chrTest1", 0, 1, '+', genome) == 'A'

    # Should get the same result if we use 1-based coords
    assert tlr.fetch_seq("chrTest1", 1, 1, '+', genome, indexing = 1) == 'A'

def test_firstbase_minus_strand():
    """ Fetch the first base of this genome on the + strand.
            ACTGACTGACTGAAATAAGAAACTGACTG
        The correct answer is T
    """
    genome_file = "talon_label_reads/test_inputs/toy_genome.fa"
    genome = pyfaidx.Fasta(genome_file, sequence_always_upper=True,
                           one_based_attributes=False)
    assert tlr.fetch_seq("chrTest1", 0, 1, '-', genome) == 'T'

    # Should get the same result if we use 1-based coords
    assert tlr.fetch_seq("chrTest1", 1, 1, '-', genome, indexing = 1) == 'T'

def test_range():
    """
          01   5    10   15   20   25
          ACTGACTGACTGAAATAAGAAACTGACTG (+)

          The goal here is to get sequence 'TGAAATAAGA'
          In 0-based coordinates, this is chrTest1:10-20.
    """
    genome_file = "talon_label_reads/test_inputs/toy_genome.fa"
    genome = pyfaidx.Fasta(genome_file, sequence_always_upper=True,
                           one_based_attributes=False)
    assert tlr.fetch_seq("chrTest1", 10, 20, '+', genome) == 'TGAAATAAGA'

    # Opposite strand
    assert tlr.fetch_seq("chrTest1", 10, 20, '-', genome) == 'TCTTATTTCA'

