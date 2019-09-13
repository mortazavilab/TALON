# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# create_GTF_abundance_from_database.py is designed to generate a GTF
# as well as an abundance file with the same filtering options. 

from optparse import OptionParser
import subprocess
import os
import sys
import shlex
from .create_abundance_file_from_database import main as create_abundance_file_main
from .filter_talon_transcripts import main as filter_transcripts_main
from .create_GTF_from_database import main as gtf_from_db_main

parser = OptionParser(description="""A script to generate a GTF and abundance file
							with the same filtering options.""")

parser.add_option("--db", dest = "database",
    help = "TALON database", metavar = "FILE", type = "string")
parser.add_option("--annot", "-a", dest = "annot",
    help = """Which annotation version to use. Will determine which
              annotation transcripts are considered known or novel
              relative to. Note: must be in the TALON database.""",
    type = "string")
parser.add_option("--build", "-b", dest = "build",
    help = "Genome build to use. Note: must be in the TALON database.",
    type = "string")
parser.add_option("--filter", dest ="filtering", action='store_true',
                  help = "If this option is set, the transcripts in the  \
                  database will be filtered prior to GTF creation \
                  (for more information, see filter_talon_transcripts.py)")
parser.add_option("--pairings", "-p",  dest = "pairings_file",
    help = """Optional (only relevant if filter = true): A file indicating
              which datasets should be considered together when filtering
              novel transcripts (i.e. biological replicates).
              Format: Each line of the file constitutes a group, with
              member datasets separated by commas.
              If no file is provided, then novel transcripts appearing in
              any two datasets will be accepted.""",
    metavar = "FILE", type = "string", default = None)
parser.add_option("--o", dest = "outprefix", help = "Prefix for output file",
    metavar = "FILE", type = "string")

(opt, args) = parser.parse_args()

# post-TALON_tools dir
tpath = os.path.dirname(os.path.realpath(sys.argv[0]))
print(tpath)

# make abundance file
db = opt.database
annot = opt.annot
build = opt.build 
o = opt.outprefix
filtering = opt.filtering
pairings = opt.pairings_file
create_abundance_file_arguments = (
    "--db {} --annot {} --build {} --o {}".format(tpath, db, annot, build, o))

if filtering:
	create_abundance_file_arguments+= ' --filter'
	if pairings != None:
		create_abundance_file_arguments+= ' --pairings {}'.format(pairings)

# TODO: Call a function with argument instead of using argv and calling main
sys.argv = ["create_abundance_file_from_database.py"] + shlex.split(
	create_abundance_file_arguments)
create_abundance_file_main()

# make whitelist file for GTF
if filtering:
	outfile = o+'_whitelist'
	filter_arguments = "--db {} --annot {} --o {}".format(
		tpath, db, annot, outfile)
	if pairings != None:
		filter_arguments+=' --pairings {}'.format(pairings)
	sys.argv = ["filter_talon_transcripts.py"] + shlex.split(filter_arguments)
	filter_transcripts_main()

# make GTF
gtf_from_db_arguments = "--db {} --build {} --annot {} --o {}".format(
	tpath, db, build, annot, o)
if filtering:
	gtf_from_db_arguments +=' --whitelist {}'.format(outfile)
else:
	gtf_from_db_arguments +=' --observed'
	# pfile = open(pairings, 'r')
	# pairing_str = pfile.read()
	# pfile.close()
	# pairing_str.replace(',', '\n')
	# ofile = o+'_datasets'
	# ofile = open(ofile, 'w')
	# ofile.write(pairing_str)
	# cmd+=' --datasets {}'.format(o+'_datasets')
sys.argv = ["create_GTF_from_database.py"] + shlex.split(gtf_from_db_arguments)
gtf_from_db_main()
