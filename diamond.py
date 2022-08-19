import shlex
import subprocess
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Tuple, Set

import pandas as pd

from taxonomy import Taxonomy
from utils import FilePathType


def run_diamond_blastp(
    cleaned_fastafile: FilePathType,
    diamond_output: FilePathType,
    database: FilePathType,
    threads: int,
) -> None:
    """Run diamond blastP

    Args:
        cleaned_fastafile (FilePathType): fasta database with headers that have no spaces
        diamond_output (FilePathType): name of output diamond file
        database (FilePathType): diamond database file
        threads (int): number of worker threads
    """
    user_args = (
        f"-d {database} -q {cleaned_fastafile} -o {diamond_output} --threads {threads}"
    )
    common_args = f"--max-target-seqs 3 -e 1e-5 --query-cover 10 --subject-cover 10 --quiet -f 6 qseqid sseqid pident evalue bitscore"
    command = f"diamond blastp {user_args} {common_args}"
    subprocess.run(shlex.split(command))


def preprocess_diamond_output(
    diamond_output: FilePathType, new_output: FilePathType
) -> None:
    """Diamond automatically truncates fasta headers after the first space,
    but that was accounted for in a preprocessing step to clean the fasta
    headers. This undoes that temporary name change and keeps this output
    as the de facto diamond output.

    Args:
        diamond_output (FilePathType): name of original diamond output
        new_output (FilePathType): name of new diamond output that has
            the query names unmangled
    """
    with open(diamond_output) as infile, open(new_output, "w") as outfile:
        for line in infile:
            outfile.write(line.replace("~~", " "))


def parse_diamond_hits(
    diamond_output: FilePathType,
    taxonomy: Dict[str, Taxonomy],
    pident_cutoff: float = 0.2,
    protein_level: bool = False,
) -> DefaultDict[str, List[Taxonomy]]:
    """Parse the diamond output file. First, the subject hit accession number is mapped
    to a `Taxonomy` object that was read from a file. Then, the query genome is mapped
    to a list of all `Taxonomy` hits.

    Args:
        diamond_output (FilePathType): name of cleaned diamond output
        taxonomy (Dict[str, Taxonomy]): taxonomy mapping from genome accessions to `Taxonomy`
            objects --> read from a provided file in the main script
        pident_cutoff (float, optional): protein identity cutoff for counting a hit. 
            Defaults to 0.2.
        protein_level (bool, optional): get taxonomy at protein level. 
            Defaults to False.

    Returns:
        DefaultDict[str, List[Taxonomy]]: maps query genome to a list of 
            all `Taxonomy` hits
    """
    diamond_hits: DefaultDict[str, List[Taxonomy]] = defaultdict(list)
    with open(diamond_output) as infile:
        for line in infile:
            # ignore evalue and bitscore for now
            ptn_query, ptn_subject, _pident, *_ = line.rstrip().split("\t")
            pident = float(_pident)
            if pident >= pident_cutoff:
                reference = ptn_subject.rsplit("_", 1)[0]
                ref_tax = taxonomy[reference]
                # get taxonomy of genome or protein
                query = ptn_query if protein_level else ptn_query.rsplit("_", 1)[0]
                diamond_hits[query].append(ref_tax)
    return diamond_hits


def taxonomy_hits_to_df(diamond_hits: DefaultDict[str, List[Taxonomy]]) -> pd.DataFrame:
    """Convert a dictionary of a list of `Taxonomy` hits into a pandas dataframe.
    All downstream operations are best performed with dataframe logic, so this is 
    helping intermediate step.

    Args:
        diamond_hits (DefaultDict[str, List[Taxonomy]]): maps query genome to a list of 
            all `Taxonomy` hits

    Returns:
        pd.DataFrame: converted diamond_hits dictionary into a dataframe
    """
    _df: List[Dict[str, str]] = list()
    for query, hits in diamond_hits.items():
        for hit in hits:
            _df.append({"Genome": query, **asdict(hit)})
    return pd.DataFrame(_df)


def top_taxa(
    diamond_hits: pd.DataFrame,
    order: bool = True,
    family: bool = True,
    subfamily: bool = True,
    genus: bool = True,
) -> pd.DataFrame:
    """Calculate the number of taxa representatives at the included
    taxonomic levels. A representative is an entire taxonomic classification
    grouping. Ie "Caudovirales, Myoviridae, Tevenvirinae" is a representative.

    Args:
        diamond_hits (pd.DataFrame): dataframe of taxonomy hits for all proteins
            from the query genomes
        order (bool, optional): include order. Defaults to True.
        family (bool, optional): include family. Defaults to True.
        subfamily (bool, optional): include subfamily. Defaults to True.
        genus (bool, optional): include genus. Defaults to True.

    Returns:
        pd.DataFrame: counts of all unique taxa representatives for each query
            genome
    """
    base_groups = ["Genome", "order", "family", "subfamily", "genus"]
    masks = [True, order, family, subfamily, genus]
    chosen_groups = [group for group, mask in zip(base_groups, masks) if mask]
    return diamond_hits.groupby(by=chosen_groups, as_index=False).size()


def top_individual_taxa(diamond_hits: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Counterpart function to `top_taxa`. Calculate the number of overall
    hits for each taxonomic level (excluding genus), using the `top_taxa` function.

    This does not consider entire taxonomic representatives. For example:

    Caudovirales    Myoviridae  Tevenvirinae    10
    Caudovirales    Myoviridae  Peduovirinae    2

    Will yield results for the given genome:
    Caudovirales    12
    Myoviridae      12
    Tevenvirinae    10
    Peduovirinae    2

    Args:
        diamond_hits (pd.DataFrame): dataframe of taxonomy hits for all proteins
            from the query genomes

    Returns:
        Dict[str, pd.DataFrame]: maps taxonomic level to the respective
            top_taxa dataframe that only considered that level
    """
    levels = ("order", "family", "subfamily", "genus")

    top_indiv_taxa: Dict[str, pd.DataFrame] = dict()
    for taxon in ("order", "family", "subfamily"):
        taxon_mask = {level: level == taxon for level in levels}
        top = (
            top_taxa(diamond_hits, **taxon_mask)
            .assign(
                max_size=lambda df: df.groupby("Genome", as_index=False)[
                    "size"
                ].transform(max)
            )
            .query("size == max_size")
            .drop("max_size", axis=1)
        )
        top_indiv_taxa[taxon] = top

    return top_indiv_taxa


def top_2_taxa(
    genome: str,
    level: str,
    tophits: pd.DataFrame,
    genome_topreps: pd.DataFrame,
    current_pred: str,
) -> Set[str]:
    """Top 2 taxa predictions are the second representative
    or the top overall prediction without considering an 
    entire representative

    Args:
        genome (str): query genome
        level (str): taxonomic level: {'order', 'family', 'subfamily'}
        tophits (pd.DataFrame): top taxonomic hits for the above level
        # TODO: should change tophits to be like genome_topreps....
        genome_topreps (pd.DataFrame): top taxa representatives for the given genome
        current_pred (str): current taxon prediction

    Returns:
        Set[str]: set of at most 2 possible choices of top taxa, which will
            be compared to the current prediction
    """
    level_columns = {"order": 1, "family": 2, "subfamily": 3}
    col = level_columns[level]

    # get taxa of second top representative
    second_size: int = int(genome_topreps.iloc[1, 4])
    second_reps = genome_topreps.query("size == @second_size")
    if len(second_reps) == 1:
        # no ties for second predicted taxa
        second_toprep: str = str(genome_topreps.iloc[1, col])
    else:
        # since there are ties for a second predicted taxa,
        # this is likely an ambiguous case,
        # so choose the other genome
        second_reps_filtered = second_reps.loc[second_reps[level] != current_pred]
        second_toprep: str = str(
            second_reps_filtered.iloc[0, col]
        ) if not second_reps_filtered.empty else "ambiguous"

    # get the overall top taxa regardless of true taxa combinations seen with reps
    genome_tophits = tophits.query("Genome == @genome")
    if len(genome_tophits) == 1:
        # only 1 tophit at this level means just choose that taxa
        overall_top: str = tophits.set_index("Genome").at[genome, level]
    else:
        # if there are more than 2 tophits for a given taxa, then choose
        # the one that is consistent with second_toprep
        for tophit in genome_tophits[level]:
            if tophit == second_toprep:
                overall_top = tophit
                break
        else:
            overall_top = "ambiguous"

    return {overall_top, second_toprep}


def make_taxa_ambiguous(genome: str, topreps: pd.DataFrame, level: str) -> None:
    """Update the taxon for a given genome to `ambiguous`

    Args:
        genome (str): query genome
        topreps (pd.DataFrame): current top taxa representatives dataframe
        level (str): taxonomic level: {'order', 'family', 'subfamily'}
    """
    topreps.loc[topreps["Genome"] == genome, level] = "ambiguous"


def profile_taxa(diamond_hits: pd.DataFrame, unique_species: Set[str]) -> pd.DataFrame:
    """_summary_

    Args:
        diamond_hits (pd.DataFrame): _description_
        unique_species (Set[str]): _description_

    Returns:
        pd.DataFrame: _description_
    """
    topreps = top_taxa(diamond_hits, genus=False)
    top_indiv_taxa = top_individual_taxa(diamond_hits)

    result = (
        topreps.sort_values(by="size", ascending=False)
        .drop_duplicates(subset="Genome", keep="first")
        .reset_index(drop=True)
    )

    for row in result.itertuples():
        genome_topreps = topreps.query("Genome == @row.Genome").sort_values(
            by="size", ascending=False
        )
        # the taxon will stay the same if there is only 1 reported taxa for this phage
        n_orders, n_families, n_subfamilies = [
            genome_topreps[level].unique().size
            for level in ("order", "family", "subfamily")
        ]
        if not n_orders == 1:
            ######## ORDER ##########
            # means there are more to choose from, meaning that the order
            # may not be the overall_top_order

            top_orders = top_2_taxa(
                row.Genome, "order", top_indiv_taxa["order"], genome_topreps, row.order
            )

            # since there are multiple orders to choose from,
            # if the indicated order is either the top overall order
            # or the same as the second next order, then keep it the same
            # otherwise, the taxonomy is too ambiguous
            if row.order not in top_orders:
                # all lower taxa must be ambiguous too
                for level in ("order", "family", "subfamily"):
                    make_taxa_ambiguous(row.Genome, result, level)
                continue

        ####### FAMILY #######
        # same logic as above
        if not n_families == 1:
            top_families = top_2_taxa(
                row.Genome,
                "family",
                top_indiv_taxa["family"],
                genome_topreps,
                row.family,
            )
            if row.family not in top_families:
                for level in ("family", "subfamily"):
                    make_taxa_ambiguous(row.Genome, result, level)
                continue

        ####### SUBFAMILY #######
        if not n_subfamilies == 1:
            top_subfamilies = top_2_taxa(
                row.Genome,
                "subfamily",
                top_indiv_taxa["subfamily"],
                genome_topreps,
                row.subfamily,
            )
            if row.subfamily not in top_subfamilies:
                make_taxa_ambiguous(row.Genome, result, "subfamily")

    result.loc[result["size"] < 3, ["order", "family", "subfamily"]] = "unknown"
    result = result.drop("size", axis=1)

    included_species = set(result["Genome"].unique())

    AMBIGUOUS_TAXA = {level: "ambiguous" for level in ("order", "family", "subfamily")}

    excluded_species: List[Dict[str, str]] = list()
    for species in unique_species:
        if species not in included_species:
            row = {"Genome": species, **AMBIGUOUS_TAXA}
            excluded_species.append(row)
    return pd.concat([result, pd.DataFrame(excluded_species)])


def most_signif_hits(diamond_output: FilePathType) -> pd.DataFrame:
    """Filter the raw diamond output to choose the hit that is 
    most significant for each query. Most significant is defined as the lowest evalue.

    Args:
        diamond_output (str): path to raw diamond output

    Returns:
        pd.DataFrame: dataframe of most significant hits for every query protein
    """
    return (
        pd.read_table(
            diamond_output, names=["ptn_query", "subject", "pident", "evalue", "score"]
        )
        .sort_values(by="evalue", ascending=True)
        .drop_duplicates(subset="ptn_query", keep="first")
        .assign(subject_genome=lambda df: df["subject"].str.rsplit("_", 1).str[0])
        .sort_values(by="ptn_query")
        .reset_index(drop=True)
    )


def best_ptn_taxonomy(df: pd.DataFrame, taxonomy: Dict[str, Taxonomy]) -> pd.DataFrame:
    """Helper funcciton to get the most significant taxonomy for every query protein.
    Converts column of `Taxonomy` objects into separate columns.

    Args:
        df (pd.DataFrame): dataframe of most significant hits for each query protein.
            Calculated from `most_signif_hits`.
        taxonomy (Dict[str, Taxonomy]): taxonomy mapping from genome accessions to `Taxonomy`
            objects --> read from a provided file in the main script

    Returns:
        pd.DataFrame: most significant taxa for each query protein
    """

    def _taxonomy_to_dict(genome: str, taxonomy: Dict[str, Taxonomy]):
        return taxonomy[genome].to_dict()

    # unwrap `Taxonomy` objects into individual columns
    _taxa = pd.DataFrame(
        df["subject_genome"].apply(_taxonomy_to_dict, taxonomy=taxonomy).values.tolist()
    )

    # add unwrapped `Taxonomy` columns back onto the original dataframe
    taxa = pd.concat([df, _taxa], axis=1).drop(
        ["subject", "pident", "evalue", "score", "subject_genome"], axis=1
    )
    return taxa


def diamond_main(
    _diamond_output: Path,
    cleaned_fastafile: Path,
    database: Path,
    taxonomy: Dict[str, Taxonomy],
    threads: int,
    unique_species: Set[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    run_diamond_blastp(cleaned_fastafile, _diamond_output, database, threads)
    diamond_output = Path(f"{_diamond_output}.tsv")
    preprocess_diamond_output(_diamond_output, diamond_output)

    ptn_taxonomy = most_signif_hits(diamond_output).pipe(
        best_ptn_taxonomy, taxonomy=taxonomy
    )
    diamond_hits = taxonomy_hits_to_df(parse_diamond_hits(diamond_output, taxonomy))

    # TODO: need to check query genomes that had no homologs
    # TODO: nede to add ambiguous columns to those....
    # TODO: this shouldn't be part of the of this step? during subsampling?? idk
    viral_taxonomy = profile_taxa(diamond_hits, unique_species)

    return ptn_taxonomy, viral_taxonomy
