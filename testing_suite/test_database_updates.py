import pytest
import os
from talon import talon
from talon import init_refs
from .helper_fns import get_db_cursor
@pytest.mark.dbunit

class TestDatabaseUpdates(object):

    def test_abundance(self):
        """ Try to add abundance entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        abundance = [ ( 1, "test", 5),
                      ( 2, "test", 1),
                      ( 3, "test", 2)]
        batch_size = 2
        talon.batch_add_abundance(cursor, abundance, batch_size)

        # Test if items are there
        query = "SELECT * FROM abundance"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 3
        conn.close()

    def test_datasets(self):
        """ Try to add dataset metadata to database """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        datasets = [ ( 1, "toy", "toy", "toy") ]
        talon.add_datasets(cursor, datasets)

        # Test if items are there
        query = "SELECT * FROM dataset"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 1
        conn.close()

    def test_observed(self):
        """ Try to add observed entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        observed = [ ( 1, 1, 1, "read1", "test", 1, 2, 1, 1, 0, 0, 100),
                     ( 2, 1, 1, "read2", "test", 1, 2, 1, 1, 0, 0, 100),
                     ( 3, 1, 1, "read3", "test", 1, 2, 1, 1, 0, 0, 100),
                     ( 4, 1, 8, "read4", "test", 35, 36, 32, 32,  None, None, 100) ]

        # Write observed to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/observed.tsv", 'w') as f:
            for obs in observed:
                f.write("\t".join([str(x) for x in obs]) + "\n")
                

        batch_size = 1
        talon.batch_add_observed(cursor, "scratch/db_updates/observed.tsv", batch_size)

        # Test if items are there
        query = "SELECT * FROM observed"
        cursor.execute(query)
        results = cursor.fetchall()
        assert len(results) == 4

        # Test that the 'None' values are properly recorded
        for transcript in results:
            if transcript['read_name'] == "read4":
                assert transcript['start_delta'] == None
                assert transcript['end_delta'] == None
        conn.close()

    def test_gene_annot(self):
        """ Try to add gene annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]

        # Write observed to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/gene_annot.tsv", 'w') as f:
            for entry in annot:
                f.write("\t".join([str(x) for x in entry]) + "\n")

        batch_size = 1
        talon.batch_add_annotations(cursor, "scratch/db_updates/gene_annot.tsv", 
                                    "gene", batch_size)

        # Test if items are there
        query = "SELECT * FROM gene_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_transcript_annot(self):
        """ Try to add transcript annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/transcript_annot.tsv", 'w') as f:
            for entry in annot:
                f.write("\t".join([str(x) for x in entry]) + "\n")

        batch_size = 2
        talon.batch_add_annotations(cursor, "scratch/db_updates/transcript_annot.tsv", 
                                    "transcript", batch_size)

        # Test if items are there
        query = "SELECT * FROM transcript_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_exon_annot(self):
        """ Try to add exon annotation entries to database in batches
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"

        annot = [ ( 1, "toy", "TALON", "status", "NOVEL"),
                  ( 2, "toy", "TALON", "status", "NOVEL") ]

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/exon_annot.tsv", 'w') as f:
            for entry in annot:
                f.write("\t".join([str(x) for x in entry]) + "\n")

        batch_size = 3
        talon.batch_add_annotations(cursor, "scratch/db_updates/exon_annot.tsv", 
                                    "exon", batch_size)

        # Test if items are there
        query = "SELECT * FROM exon_annotations WHERE value = 'NOVEL'"
        cursor.execute(query)
        assert len(cursor.fetchall()) == 2
        conn.close()

    def test_gene_update(self):
        """ Try to add novel gene entries to database while ignoring duplicates
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        run_info = talon.init_run_info(database, build)
        talon.get_counters(database)

        init_refs.make_temp_novel_gene_table(cursor, build)
        talon.create_gene("chr4", 1, 1000, "+", cursor, "temp_gene")

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/genes.tsv", 'w') as f:
            cursor.execute("SELECT gene_ID, strand FROM temp_gene")
            for entry in cursor.fetchall():
                f.write("\t".join([str(x) for x in entry]) + "\n")

        talon.batch_add_genes(cursor, "scratch/db_updates/genes.tsv", 10)

        # Test if gene with ID 6 is there, but make sure we didn't add 
        # duplicates of the other genes
        query = "SELECT * FROM genes"
        gene_IDs = [ x['gene_ID'] for x in cursor.execute(query)]
        assert 7 in gene_IDs
        assert len(gene_IDs) == 7
        conn.close()

    def test_transcript_update(self):
        """ Try to add novel transcript entries to database while ignoring 
            duplicates
        """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        transcript_dict = init_refs.make_transcript_dict(cursor, build)
        database = "scratch/toy.db"
        talon.get_counters(database)
        talon.create_transcript("chr1", 1, 1000, 1, (1,), (1,2), transcript_dict)

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/transcripts.tsv", 'w') as f:
            for transcript in transcript_dict.values():
                if type(transcript) is dict:
                    entry = "\t".join([ str(x) for x in ( transcript['transcript_ID'],
                                                          transcript['gene_ID'],
                                                          transcript['start_exon'],
                                                          transcript['jn_path'],
                                                          transcript['end_exon'],
                                                          transcript['start_vertex'],
                                                          transcript['end_vertex'],
                                                          transcript['n_exons'] ) ])
                    f.write(entry + "\n")

        batch_size = 5
        talon.batch_add_transcripts(cursor, "scratch/db_updates/transcripts.tsv", batch_size)

        # Test if transcript with ID 8 is there, but make sure we didn't add
        # duplicates of the others
        query = "SELECT * FROM transcripts"
        cursor.execute(query)
        transcripts = cursor.fetchall()
        transcript_IDs = [ x['transcript_ID'] for x in transcripts]
        assert 8 in transcript_IDs
        assert len(transcript_IDs) == 8

        # Test if None value was handled correctly
        for transcript in transcripts:
            if transcript['transcript_ID'] == 8:
                assert transcript['jn_path'] == None

        conn.close()

    def test_edge_update(self):
        """ Try to add novel exons and introns. """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        edge_dict = init_refs.make_edge_dict(cursor)
        run_info = talon.init_run_info(database, build)
        orig_n_edges = talon.edge_counter.value()

        talon.create_edge(2, 1, "exon", "-", edge_dict)

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/edges.tsv", 'w') as f:
            for edge in list(edge_dict.values()):
                if type(edge) is dict:
                    entry = "\t".join([str(x) for x in [edge['edge_ID'], edge['v1'],
                                                        edge['v2'], edge['edge_type'],
                                                        edge['strand']]] )
                    f.write(entry + "\n")

        batch_size = 10
        talon.batch_add_edges(cursor, "scratch/db_updates/edges.tsv", batch_size)
        
        # Test if the edge table has the correct number of edges now
        query = "SELECT * FROM edge"
        cursor.execute(query)
        edge_IDs = [ x['edge_ID'] for x in cursor.fetchall()]
        assert orig_n_edges + 1 in edge_IDs
        assert len(edge_IDs) == orig_n_edges + 1
        conn.close()

    def test_location_update(self):
        """ Update locations """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"
        talon.get_counters(database)
        location_dict = init_refs.make_location_dict(build, cursor)
        run_info = talon.init_run_info(database, build)
        orig_n_pos = talon.vertex_counter.value()

        talon.create_vertex("chr4", 2000, location_dict, run_info)
   
        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/loc.tsv", 'w') as f:
            for chrom_dict in location_dict.values():
                for loc in list(chrom_dict.values()):
                    if type(loc) is dict:
                        entry = ("\t".join([ str(x) for x in (loc['location_ID'],
                                                            loc['genome_build'],
                                                            loc['chromosome'],
                                                            loc['position'])]))
                        f.write(entry + "\n")

        batch_size = 10
        talon.batch_add_locations(cursor,"scratch/db_updates/loc.tsv" , batch_size)

        # Test if the table has the correct number of locations now
        query = "SELECT * FROM location"
        cursor.execute(query)
        loc_IDs = [ x['location_ID'] for x in cursor.fetchall()]
        assert orig_n_pos + 1 in loc_IDs
        assert len(loc_IDs) == orig_n_pos + 1
        conn.close()

    def test_vertex2gene_update(self):
        """ Update vertex to gene relationships """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        vertex_2_gene = init_refs.make_vertex_2_gene_dict(cursor)

        talon.update_vertex_2_gene(2, (1,2), "-", vertex_2_gene) 
        talon.update_vertex_2_gene(1, (1,2,3,4,5,6), "+", vertex_2_gene)

        # Write to file
        os.system("mkdir -p scratch/db_updates/")
        with open("scratch/db_updates/v2g.tsv", 'w') as f: 
            for vertex_ID, gene_set in vertex_2_gene.items():
                for gene in gene_set:
                    entry = "\t".join([ str(x) for x in (vertex_ID, gene[0])])
                    f.write(entry + "\n")

        batch_size = 100
        talon.batch_add_vertex2gene(cursor, "scratch/db_updates/v2g.tsv", batch_size)

        # Use queries to check if the insert worked as expected
        query = "SELECT * FROM vertex WHERE vertex_ID = '1'"
        cursor.execute(query)
        gene_IDs = [ x['gene_ID'] for x in cursor.fetchall()]
        assert gene_IDs == [1, 2]

        query = "SELECT * FROM vertex WHERE gene_ID = '1'"
        cursor.execute(query)
        vertex_IDs = [ x['vertex_ID'] for x in cursor.fetchall()]
        assert vertex_IDs == [1, 2, 3, 4, 5, 6]

    def test_counter_update(self):
        """ Update counters """
        conn, cursor = get_db_cursor()
        build = "toy_build"
        database = "scratch/toy.db"

        talon.get_counters(database)

        # Change the counter values to some arbitrary numbers
        for i in range(10): talon.gene_counter.increment()
        for i in range(20): talon.transcript_counter.increment()
        for i in range(2): talon.edge_counter.increment()
        for i in range(5): talon.vertex_counter.increment()
        n_datasets = 30
        for i in range(6): talon.observed_counter.increment()

        # Now try the update
        talon.update_counter(cursor, n_datasets)

        # Check results with queries
        cursor.execute("""SELECT * FROM counters""")
        for category,value in cursor.fetchall():
            if category == "genes":
                assert value == 16
            elif category == "transcripts":
                assert value == 27
            elif category == "edge":
                assert value == 33
            elif category == "vertex":
                assert value == 39
            elif category == "observed":
                assert value == 6
            elif category == "dataset":
                assert value == 30
            else:
                if category != "genome_build":
                    pytest.fail("Unexpected entry in counters table")

        conn.close()
