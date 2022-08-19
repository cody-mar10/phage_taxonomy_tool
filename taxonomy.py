from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict


@dataclass(eq=True, frozen=True)
class Taxonomy:
    """Stores order, family, subfamily, and genus for viral taxonomy"""

    order: str
    family: str
    subfamily: str
    genus: str

    def to_dict(self) -> Dict[str, str]:
        return asdict(self)


def parse_taxonomy(taxonomy: Path) -> Dict[str, Taxonomy]:
    """Parse the taxonomy file into a dictionary that maps
    the genome name to a `Taxonomy` dataclass.

    Args:
        taxonomy (Path): path to taxonomy file

    Returns:
        Dict[str, Taxonomy]: maps the genome name to a `Taxonomy` 
            dataclass
    """
    with taxonomy.open() as fp:
        # skip header
        fp.readline()
        taxa: Dict[str, Taxonomy] = dict()
        for line in fp:
            accession, *taxons = line.rstrip().split("\t")
            taxa[accession] = Taxonomy(*taxons)
        return taxa
