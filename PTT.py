#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2020

##############
## This script is associated with the manuscript
## Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages
## Kieft and Zhou et al. 2020
##############

# Usage: $ python3 PTT.py -i <input_fasta_file> -t <threads> -f <format>
# PTT: phage taxonomy tool
# Version comment: the database was compiled using Diamond v2.0.0.138

# Contact Kristopher Kieft (kieft@wisc.edu) with questions regarding this script

"""
2022-08-15 Cody Martin - i am heavily modifying this script to work better
so now the version is 0.1.1
"""

import argparse
import datetime
from pathlib import Path
from typing import Set

from __init__ import __version__
from diamond import diamond_main
from taxonomy import parse_taxonomy
from utils import FilePathType, log_command, preprocess_ptn_fasta, remove_tmp_files


def parse_args() -> argparse.Namespace:
    PTT_path = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Estimates phage taxonomy.")
    parser.add_argument("--version", action="version", version=f"PTT v{__version__}")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help='input protein fasta file. Must be in format "name_#" with no spaces in protein name (Prodigal format).',
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="number of threads [default = %(default)s]",
    )
    parser.add_argument(
        "-f",
        "--format",
        default="prot",
        choices={"prot", "nucl"},
        help="format of input [default = %(default)s]",
    )
    parser.add_argument(
        "-d",
        "--database",
        default=PTT_path.joinpath("PTT_database.dmnd").as_posix(),
        help="path to diamond database (default: %(default)s)",
    )
    parser.add_argument(
        "-x",
        "--taxonomy",
        default=PTT_path.joinpath("PTT_virus_taxonomy.tsv").as_posix(),
        help="path to reference taxonomy file (default: %(default)s)",
    )
    return parser.parse_args()


# subprocess.run("rm " + str(base) + ".PTT.diamond.out 2>/dev/null", shell=True)
# subprocess.run("rm " + str(base) + ".temp.fasta 2>/dev/null", shell=True)
# subprocess.run("rm " + str(base) + ".temp.faa 2>/dev/null", shell=True)

# TODO: move decorator to main so it can take in the output directory name
def PTT_main(
    fastafile: FilePathType,
    database: FilePathType,
    taxonomyfile: FilePathType,
    threads: int,
    format: str,
    today: datetime.date,
    outdir: Path,
):
    fastafile = Path(fastafile)
    base_fastafile = fastafile.stem
    database = Path(database)
    taxonomyfile = Path(taxonomyfile)
    PTT_output_base = outdir.joinpath(f"{base_fastafile}.PTT")

    tax_map = parse_taxonomy(taxonomyfile)

    tmpfiles: Set[FilePathType] = set()

    cleaned_fastafile = outdir.joinpath(f"{base_fastafile}.tmp.faa")
    tmpfiles.add(cleaned_fastafile)
    unique_sp = preprocess_ptn_fasta(fastafile, cleaned_fastafile)

    _diamond_output = Path(f"{PTT_output_base}.diamond")
    tmpfiles.add(_diamond_output)
    ptn_taxonomy, viral_taxonomy = diamond_main(
        _diamond_output, cleaned_fastafile, database, tax_map, threads, unique_sp
    )

    ptn_taxonomy.to_csv(
        f"{PTT_output_base}.protein-taxonomy.tsv", sep="\t", index=False
    )
    viral_taxonomy.to_csv(
        f"{PTT_output_base}.viral-taxonomy.tsv", sep="\t", index=False
    )

    remove_tmp_files(*tmpfiles)


if __name__ == "__main__":
    args = parse_args()
    fastafile = Path(args.input)
    database = Path(args.database)
    taxonomyfile = Path(args.taxonomy)
    threads: int = args.threads
    format: str = args.format

    today = datetime.date.today()
    outdir = Path(f"PTT_{fastafile.stem}_{today}")
    if not outdir.exists():
        outdir.mkdir()
    logfile = outdir.joinpath("PTT.log").as_posix()
    log_command(logfile)(PTT_main)(
        fastafile, database, taxonomyfile, threads, format, today, outdir
    )
