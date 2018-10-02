import pytest
import sys
sys.path.append("..")
from transcript import *
from edge import *

class TestAddEdgeToTranscript(object):
    def test1(self):
        """
        Task: 
        Add three exons to the transcript that belong in the following order:
        e1, e2, e3
        (10,20) (30,40) (45,50)
        In this version of the test, we attempt to add them in order e1, e2, e3
        """
        start = 10
        end = 50
        e1_start = 10
        e1_end = 20
        e2_start = 30
        e2_end = 40
        e3_start = 45
        e3_end = 50

        t = Transcript("t1", "chr1", start, end, "+", "gene1", None)
        e1 = Edge("e1", "chr1", e1_start, e1_end, "+", "gene1", "t1", None)
        e2 = Edge("e2", "chr1", e2_start, e2_end, "+", "gene1", "t1", None) 
        e3 = Edge("e3", "chr1", e3_start, e3_end, "+", "gene1", "t1", None)

        t.add_exon(e1)
        t.add_exon(e2)
        t.add_exon(e3)

        ids_in_order = [exon.identifier for exon in t.exons]
        assert ids_in_order == ["e1", "e2", "e3"]

    def test2(self):
        """
        Task:
        Add three exons to the transcript that belong in the following order:
        e1, e2, e3
        (10,20) (30,40) (45,50)
        In this version of the test, we attempt to add them in order e3, e2, e1
        """
        start = 10
        end = 50
        e1_start = 10
        e1_end = 20
        e2_start = 30
        e2_end = 40
        e3_start = 45
        e3_end = 50

        t = Transcript("t1", "chr1", start, end, "+", "gene1", None)
        e1 = Edge("e1", "chr1", e1_start, e1_end, "+", "gene1", "t1", None)
        e2 = Edge("e2", "chr1", e2_start, e2_end, "+", "gene1", "t1", None)
        e3 = Edge("e3", "chr1", e3_start, e3_end, "+", "gene1", "t1", None)

        t.add_exon(e3)
        t.add_exon(e2)
        t.add_exon(e1)

        ids_in_order = [exon.identifier for exon in t.exons]
        assert ids_in_order == ["e1", "e2", "e3"]

    def test3(self):
        """
        Task:
        Add three exons to the transcript that belong in the following order:
        e1, e2, e3
        (10,20) (30,40) (45,50)
        In this version of the test, we attempt to add them in order e2, e1, e3
        """
        start = 10
        end = 50
        e1_start = 10
        e1_end = 20
        e2_start = 30
        e2_end = 40
        e3_start = 45
        e3_end = 50

        t = Transcript("t1", "chr1", start, end, "+", "gene1", None)
        e1 = Edge("e1", "chr1", e1_start, e1_end, "+", "gene1", "t1", None)
        e2 = Edge("e2", "chr1", e2_start, e2_end, "+", "gene1", "t1", None)
        e3 = Edge("e3", "chr1", e3_start, e3_end, "+", "gene1", "t1", None)

        t.add_exon(e2)
        t.add_exon(e1)
        t.add_exon(e3)

        ids_in_order = [exon.identifier for exon in t.exons]
        assert ids_in_order == ["e1", "e2", "e3"]
