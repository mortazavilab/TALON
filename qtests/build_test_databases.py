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
        "--a",  "gencodeM7",
        "--l", "0",
        "--g",  "mm10", "--o", "scratch/Canx"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on Canx annotation")
