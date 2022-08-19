#!/usr/bin/env python3
from pathlib import Path
from typing import Iterator, List, Tuple, Union
import textwrap

FASTASequence = Tuple[str, str]


def fastaparser(file: Union[str, Path]) -> Iterator[FASTASequence]:
    """Parse a FASTA file by yielding (header, sequence) tuples.

    Args:
        file (str | Path): path to a valid FASTA file

    Yields:
        FASTASequence: (header, sequence) tuple
    """
    with open(file) as f:
        header = f.readline().rstrip()[1:]
        seq: List[str] = list()
        for line in f:
            line = line.rstrip()
            if line[0] == ">":
                yield header, "".join(seq)
                header = line[1:]
                seq = list()
            else:
                seq.append(line)
        yield header, "".join(seq)


def wrap_fasta(sequence: str, width: int = 75) -> Iterator[str]:
    """Wrap a sequence when outputting to a FASTA file.

    Args:
        sequence (str): biological sequence
        width (int, optional): width of a sequence line. 
            Defaults to 75 characters.

    Yields:
        str: a single sequence line
    """
    yield from textwrap.wrap(sequence, width=width)
