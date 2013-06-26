#!/usr/bin/env python

import sys
from Bio.SeqIO import parse, write
from Bio.Seq import translate, Seq
from backalign import codons


def transl(seq):
    out_seq = ""
    for codon in codons(seq):
        if codon == '---':
            aa = '-'
        elif codon == '...':
            aa = '.'
        else:
            aa = translate(codon)
        out_seq += aa
    return out_seq


for rec in parse(sys.stdin, 'fasta'):
    rec.seq = Seq(transl(rec.seq))
    write(rec, sys.stdout, 'fasta')
