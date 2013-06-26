#!/usr/bin/env python3
"""Align nucleotides to match the amino-acid alignment."""

from Bio.SeqIO import parse, write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def codons(sequence):
    """Return groups of 3 items (nucleotides) from a sequence."""
    codon = ""
    for nucl in sequence:
        codon += nucl
        if len(codon) == 3:
            yield codon
            codon = ""

def backalign(aligned_aa, unaligned_nucl):
    afn_str = ""
    codon_iter = codons(unaligned_nucl)
    for amino in aligned_aa:
        if amino == ".":
            continue
        elif amino.islower():
            next(codon_iter)
        elif amino == "-":
            afn_str += "---"
        else:
            afn_str += next(codon_iter)
    return afn_str

def main():
    for afa_rec, fn_rec in \
        zip(parse(sys.argv[1], 'fasta'), parse(sys.argv[2], 'fasta')):
        assert afa_rec.id == fn_rec.id
        afn_str = backalign(afa_rec.seq, fn_rec.seq.ungap('-'))
        write(SeqRecord(Seq(afn_str), id=afa_rec.id, description=""),
              sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
