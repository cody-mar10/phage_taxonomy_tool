import datetime
import logging
import os
import re
import sys
from pathlib import Path
from typing import Callable, Set, Union

from __init__ import __version__
from fastaparser import fastaparser, wrap_fasta

FilePathType = Union[Path, str]


def preprocess_ptn_fasta(fastafile: FilePathType, output: FilePathType) -> Set[str]:
    """Preprocess a protein fasta file by extracting only the protein name,
    which should be in the format GENOME_PTNNUM. To account for spaces that
    may be in the genome name, this temporarily replaces spaces with ~~ after
    cutting off the remaining header metadata.

    A side effect of this function is that it keeps track of all unique genomes.

    Args:
        fastafile (FilePathType): path to fasta file
        output (FilePathType): name of output processed fasta file

    Returns:
        Set[str]: set of all unique genomes
    """
    unique_species: Set[str] = set()

    # prodigal format splits everything by the ' # ' delimiter,
    # but others use a tab, and we only want the sequence name
    pattern = re.compile(" # |\t")
    with open(output, "w") as outfile:
        for header, seq in fastaparser(fastafile):
            cleaned_header = pattern.split(header, maxsplit=1)[0].replace(" ", "~~")
            unique_species.add(cleaned_header.rsplit("_", 1)[0].replace("~~", " "))
            outfile.write(f">{cleaned_header}\n")
            for seqline in wrap_fasta(seq):
                outfile.write(f"{seqline}\n")
    return unique_species


def remove_tmp_files(*files: FilePathType) -> None:
    for file in files:
        os.remove(file)


def log_command(logname: str) -> Callable:
    def wrap(func: Callable) -> Callable:
        def log(*args, **kwargs) -> None:
            logging.basicConfig(
                filename=logname, level=logging.INFO, format="%(message)s"
            )
            logging.info(f"{'Command:':10s}{' '.join(sys.argv)}")
            logging.info(f"{'Conda:':10s}{os.environ['CONDA_DEFAULT_ENV']}")
            start = datetime.datetime.now()
            logging.info(f"{'Start:':10s}{start}")

            func(*args, **kwargs)

            end = datetime.datetime.now() - start
            mins, secs = divmod(end.total_seconds(), 60)
            runtime = f"{mins} minutes and {secs:.0f} seconds"
            logging.info(f"{'Runtime:':10s}{runtime}")
            logging.info(f"{'Version:':10s}PTT v{__version__}")

            logging.info(
                "                                                               ##\n"
                "                                                             ##  ##\n"
                "                                                           ##      ##\n"
                "######   ##  ##     ##     #######   ######    #####       ##      ##\n"
                "##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##\n"
                "######   ######   ######   ##  ###   ######    ###             ##\n"
                "##       ##  ##   ##  ##   ##   ##   ##           ##           ##\n"
                "##       ##  ##   ##  ##   #######   ######   #####            ##\n"
                "                                                            #  ##  #\n"
                "                                                           # # ## # #\n"
                "                                                          #   #  #   #\n"
                "                                                         #            #\n"
            )
            logging.info("\n")

        return log

    return wrap
