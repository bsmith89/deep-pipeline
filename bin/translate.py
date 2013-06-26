#!/usr/bin/env python3
"""Translate nucleotide sequences."""

import sys
import optparse
from Bio.SeqIO import parse, write
from Bio.Seq import translate, Seq
from backalign import codons


_DEFAULT_FORMAT = 'fasta'
_DEFAULT_UNGAP = False

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


if __name__ == '__main__':
    usage = "usage: %prog [options] [NUCL_FILE]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--informat", dest="informat",
                      default=_DEFAULT_FORMAT,
                      help=("the format of the input file. "
                            "DEFAULT: {}").format(_DEFAULT_FORMAT))
    parser.add_option("-F", "--outformat", dest="outformat",
                      default=_DEFAULT_FORMAT,
                      help=("the format of the output. "
                            "DEFAULT: {}").format(_DEFAULT_FORMAT))
    parser.add_option("-u", "--ungap", dest="ungap", action="store_true",
                      default=_DEFAULT_UNGAP,
                      help=("should the input sequences have gap "
                            "characters ('-' and '.') removed before "
                            "translation? DEFAULT: {}").\
                           format(_DEFAULT_UNGAP))
    opts, args = parser.parse_args()

    if len(args) == 0:
        rec_iter = parse(sys.stdin, opts.informat)
    elif len(args) == 1:
        rec_iter = parse(args[0], opts.informat)
    else:
        parser.error("Too many arguments.")

    for rec in rec_iter:
        seq = rec.seq
        if opts.ungap:
            seq = seq.ungap('-').ungap('.')
        rec.seq = Seq(transl(seq))
        write(rec, sys.stdout, opts.outformat)
