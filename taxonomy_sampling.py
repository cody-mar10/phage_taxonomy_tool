#!/usr/bin/env python3
import argparse
import multiprocessing
from itertools import repeat
from pathlib import Path
import shutil

import pandas as pd

from diamond import parse_diamond_hits, profile_taxa, taxonomy_hits_to_df
from taxonomy import parse_taxonomy
from utils import FilePathType, log_command


def randomly_sample_proteins(
    diamond_hits: pd.DataFrame, seed: int, sample_unique_proteins: bool = False
) -> pd.DataFrame:
    n = (
        diamond_hits["Genome"].unique().size
        if sample_unique_proteins
        else len(diamond_hits)
    )
    random_sample = diamond_hits.sample(
        n=n, replace=True, random_state=seed
    ).reset_index(drop=True)
    return random_sample


def _main(
    viral_taxonomy: FilePathType,
    taxonomyfile: Path,
    diamond_output: FilePathType,
    outdir: Path,
    sample_proteins: bool,
    seed: int,
) -> None:
    with open(viral_taxonomy) as fp:
        # skip header
        fp.readline()
        unique_species = {line.split("\t", 1)[0] for line in fp}
    tax_map = parse_taxonomy(taxonomyfile)
    diamond_hits = taxonomy_hits_to_df(
        parse_diamond_hits(diamond_output, tax_map, protein_level=True)
    )
    random_sample = randomly_sample_proteins(diamond_hits, seed, sample_proteins)

    # random_sample currently has proteins in it
    # so we need to strip the protein identifier to
    # convert to the genome
    random_sample["Genome"] = random_sample["Genome"].str.rsplit("_", 1).str[0]
    results = profile_taxa(random_sample, unique_species).sort_values(by="Genome")

    output = outdir.joinpath(f"{Path(viral_taxonomy).stem}_seed{seed}.tsv")
    results.to_csv(output, sep="\t", index=False)


def main(
    viral_taxonomy: FilePathType,
    taxonomyfile: Path,
    diamond_output: FilePathType,
    outdir: Path,
    sample_proteins: bool,
    jobs: int,
    iterations: int,
) -> None:
    with multiprocessing.Pool(processes=jobs) as pool:
        pool.starmap(
            _main,
            zip(
                repeat(viral_taxonomy),
                repeat(taxonomyfile),
                repeat(diamond_output),
                repeat(outdir),
                repeat(sample_proteins),
                range(iterations),
            ),
        )


if __name__ == "__main__":
    PTT_path = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Sample diamond blastP hits for query proteins to bootstrap taxonomic profiling predictions"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=".PTT.viral-taxonomy.tsv file from the output of the initial PTT.py run on the real data",
    )
    parser.add_argument(
        "-d",
        "--diamond-output",
        required=True,
        help=".PTT.diamond.tsv file from the output of the initial PTT.py run on the real data",
    )
    parser.add_argument(
        "-b",
        "--n-bootstraps",
        type=int,
        default=10,
        help="number of boostreaps or sampling iterations to perform (default: %(default)s)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=10,
        help="number of parallel workers (default: %(default)s)",
    )
    parser.add_argument(
        "-x",
        "--taxonomy",
        default=PTT_path.joinpath("PTT_virus_taxonomy.tsv").as_posix(),
        help="path to reference taxonomy file (default: %(default)s)",
    )
    parser.add_argument(
        "--sample-proteins",
        default=False,
        action="store_true",
        help="use to sample individual proteins instead of sampling diamond hits",
    )

    args = parser.parse_args()

    results_dir = Path(args.input).parent
    outdir = results_dir.joinpath(Path("random_sampling"))

    if not outdir.exists():
        outdir.mkdir()
    else:
        shutil.rmtree(outdir)
        outdir.mkdir()

    logname = outdir.joinpath("PTT_random_sampling.log").as_posix()

    log_command(logname)(main)(
        viral_taxonomy=args.input,
        taxonomyfile=Path(args.taxonomy),
        diamond_output=args.diamond_output,
        outdir=outdir,
        sample_proteins=args.sample_proteins,
        jobs=args.jobs,
        iterations=args.n_bootstraps,
    )

