import subprocess
import sys
from pathlib import Path

R_SCRIPT_FILE = Path(__file__).parent / Path("r_scripts") / Path("generate_talon_report.R")


def main():
    args = ["Rscript", str(R_SCRIPT_FILE)] + sys.argv[1:]
    try:
        subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)
    except FileNotFoundError as e:
        print("This is a wrapper command for an R script. " "Do you have GNU R installed?\n" + str(e))


if __name__ == "__main__":
    main()
