# This script builds the databases needed for the TALON testing suite to run

import os
import subprocess
import sys
sys.path.append("filtering")
import make_minimal_db_for_filtering as mmdb

try:
   subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/toy_transcript/toy_annot.gtf",
        "--a",  "toy_annot",
        "--l", "0",
        "--g",  "toy_build", "--o", "scratch/toy"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on toy artificial annotation")

try:
   subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/toy_transcript/toy_annot.gtf",
        "--a",  "toy_annot",
        "--l", "0",
        "--g",  "toy_build", "--o", "scratch/toy_mod"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on toy artificial annotation")

try:
   subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/monoexonic_handling/monoexon.gtf",
        "--a",  "monoexon",
        "--5p", "99",
        "--3p", "99",
        "--l", "0",
        "--g",  "monoexon", "--o", "scratch/monoexon"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on monoexon artificial annotation")

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/Canx_example/Canx.gtf",
        "--a",  "gencode_vM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Canx"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Canx annotation")

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/Canx_example/Canx_refseq.gtf",
        "--a",  "gencode_vM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Canx_refseq"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Canx RefSeq annotation")

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/Map2k4_example/Map2k4.gtf",
        "--a",  "gencode_vM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Map2k4"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Map2k4 annotation")

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/chr11_and_Tcf3/chr11_and_Tcf3.gtf",
        "--a",  "gencode_vM7",
        "--5p", "500",
        "--3p", "300",
        "--idprefix", "ENCODE-mouse",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/chr11_and_Tcf3"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on chr11 annotation")

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/intergenic_GM12878/chr22.gtf",
        "--a",  "gencode_vM7",
        "--5p", "500",
        "--3p", "300",
        "--idprefix", "ENCODE-human",
        "--l", "300",
        "--g",  "hg38", "--o", "scratch/chr22"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on chr22 annotation")

try:
    os.system("mkdir -p scratch/multiexon_read_overlapping_monoexon_transcript")
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/multiexon_read_overlapping_monoexon_transcript/HMGB1P1.gtf",
        "--a",  "gencode_v29",
        "--5p", "500",
        "--3p", "300",
        "--idprefix", "ENCODEH",
        "--l", "300",
        "--g",  "hg38", "--o",
        "scratch/multiexon_read_overlapping_monoexon_transcript/talon"])

except Exception as e:
    print(e)
    sys.exit("Database initialization failed on HMGB1P1 annotation")

# Actually perform the toy_mod run
try:
    subprocess.check_output(
       ["talon",
        "--f", "input_files/toy_transcript/mod_config.csv",
        "--db", "scratch/toy_mod.db",
        "--build", "toy_build",
        "--cov", "0",
        "--identity", "0",
        "--o", "scratch/toy_mod" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on toy_mod transcripts")

# Actually perform the monoexon run
try:
    subprocess.check_output(
       ["talon",
        "--f", "input_files/monoexonic_handling/config.csv",
        "--db", "scratch/monoexon.db",
        "--build", "monoexon",
        "--cov", "0",
        "--identity", "0",
        "--o", "scratch/monoexon" ])

except Exception as e:
    print(e)
    sys.exit("TALON run failed on monoexon transcripts")

# Actually perform the chr22 TALON run
try:
    subprocess.check_output(
       ["talon",
        "--f", "input_files/intergenic_GM12878/config.csv",
        "--db", "scratch/chr22.db",
        "--build", "hg38",
        "--cov", "0",
        "--identity", "0",
        "--o", "scratch/intergenic_GM12878" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on chr22")

# Actually perform the chr11_and_Tcf3 TALON run
try:
    subprocess.check_output(
       ["talon",
        "--f", "input_files/chr11_and_Tcf3/config.csv",
        "--db", "scratch/chr11_and_Tcf3.db",
        "--build", "mm10",
        "--cov", "0",
        "--identity", "0",
        "--o", "scratch/chr11_and_Tcf3" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on chr11_and_Tcf3")

# Run the whitelist filtering script on the chr11_and_Tcf3 TALON results
try:
    subprocess.check_output(
       ["talon_filter_transcripts",
        "--db", "scratch/chr11_and_Tcf3.db",
        "-a", "gencode_vM7",
        "--datasets", "PB65_B017,PB65_B018,D12",
        "--maxFracA", "1",
        "--minCount", "1",
        "--minDatasets", "2",
        "--allowGenomic",
        "--o", "scratch/chr11_and_Tcf3_whitelist.csv" ])
except Exception as e:
    print(e)
    sys.exit("Post-TALON filtering failed on chr11_and_Tcf3")

# Run GTF script on chr11_and_Tcf3 TALON results with whitelist
try:
    subprocess.check_output(
       ["talon_create_GTF",
        "--db", "scratch/chr11_and_Tcf3.db",
        "-a", "gencode_vM7",
        "--build", "mm10",
        "--w", "scratch/chr11_and_Tcf3_whitelist.csv",
        "--o", "scratch/chr11_and_Tcf3"])
except Exception as e:
    print(e)
    sys.exit("Post-TALON GTF script failed on chr11_and_Tcf3")

# Create mock database for filtering tests
try:
    database = "scratch/filter/test.db"
    os.system("mkdir -p scratch/filter/")
    mmdb.init_mock_db(database)
except Exception as e:
    print(e)
    sys.exit("Problem creating mock database for filtering tests")


# code to get cenps-cort and rpl11-eloa only gtf
# import pyranges as pr
# df = pr.read_gtf('/Users/fairliereese/Documents/programming/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf', duplicate_attr=True).df
# gnames = ['CENPS', 'CORT', 'CENPS-CORT', 'RPL11', 'ELOA']
# df = df.loc[df.gene_name.isin(gnames)]
# df = pr.PyRanges(df)
# df.to_gtf('input_files/readthrough/readthrough.gtf')

try:
    subprocess.check_output(
       ["talon_initialize_database",
        "--f", "input_files/readthrough/readthrough.gtf",
        "--a",  "gencode_v29",
        "--5p", "500",
        "--3p", "300",
        "--idprefix", "TALON",
        "--l", "0",
        "--g",  "hg38", "--o", "scratch/readthrough"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on readthrough annotation")

# Actually perform the readthrough TALON run
try:
    subprocess.check_output(
       ["talon",
        "--f", "input_files/readthrough/config.csv",
        "--db", "scratch/readthrough.db",
        "--build", "hg38",
        "--cov", "0",
        "--identity", "0",
        "--o", "scratch/readthrough" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on readthrough")
