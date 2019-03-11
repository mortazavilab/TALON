# This script builds the databases needed for the TALON testing suite to run

import os
import subprocess
import sys

try:
   subprocess.check_output(
       ["python", "../initialize_talon_database.py",
        "--f", "input_files/toy_transcript/toy_annot.gtf",
        "--a",  "toy_annot",
        "--l", "0",
        "--g",  "toy_build", "--o", "scratch/toy"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on toy artificial annotation")

try:
    subprocess.check_output(
       ["python", "../initialize_talon_database.py",
        "--f", "input_files/Canx_example/Canx.gtf",
        "--a",  "gencode_vM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Canx"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Canx annotation")

try:
    subprocess.check_output(
       ["python", "../initialize_talon_database.py",
        "--f", "input_files/Map2k4_example/Map2k4.gtf",
        "--a",  "gencode_vM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Map2k4"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Map2k4 annotation")

try:
    subprocess.check_output(
       ["python", "../initialize_talon_database.py",
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
       ["python", "../initialize_talon_database.py",
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

# Actually perform the chr22 TALON run
try:
    subprocess.check_output(
       ["python", "../talon.py",
        "--f", "input_files/intergenic_GM12878/config.csv",
        "--db", "scratch/chr22.db",
        "--build", "hg38",
        "--o", "scratch/intergenic_GM12878" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on chr11_and_Tcf3")

# Actually perform the chr11_and_Tcf3 TALON run
try:
    subprocess.check_output(
       ["python", "../talon.py", 
        "--f", "input_files/chr11_and_Tcf3/config.csv",
        "--db", "scratch/chr11_and_Tcf3.db", 
        "--build", "mm10",
        "--o", "scratch/chr11_and_Tcf3" ])
except Exception as e:
    print(e)
    sys.exit("TALON run failed on chr11_and_Tcf3")
