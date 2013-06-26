#!/usr/bin/env python

import pandas as pd
import sys
import os
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import parse, write
import optparse



_DEFAULT_E_VALUE_CUTOFF = 1e-5
_DEFAULT_INFORMAT = 'fasta'

usage = "usage: %prog [options] HIT_FILE SEQ_FILE"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-e", "--e-value-cutoff", dest="e_value_cutoff",
                  default=_DEFAULT_E_VALUE_CUTOFF,
                  type='float',
                  help=("the maximum E-value considered a hit. "
                        "DEFAULT: {}").format(_DEFAULT_E_VALUE_CUTOFF))
parser.add_option("-f", "--informat", dest="format", default=_DEFAULT_INFORMAT,
                  help=("the format of the sequence file to be pulled "
                        "from. DEFAULT: {}").format(_DEFAULT_INFORMAT))
opts, args = parser.parse_args()

if os.path.getsize(args[0]) == 0:
    sys.stdout.write("")
    sys.exit()

hits = pd.read_table(args[0], sep='\s*', header=None)
hits = hits.select(lambda col: col in [0, 2, 4, 5, 6], axis=1)
hits.columns = ['orf', 'model', 'e_value', 'score', 'bias']
hits = hits[hits.e_value < opts.e_value_cutoff]
if len(hits) == 0:
    sys.stdout.write("")
    sys.exit()
hits['read'], hits['start'], hits['stop'] = \
        np.array(hits.orf.str.\
            match('^(.*)\(([0-9]*)-([0-9]*)\)$').values.tolist()).T
read_set = set(hits.read)
out_recs = []
for rec in parse(args[1], opts.format):
    if rec.name in read_set:
        read_hits = hits[hits.read == rec.name]
        read_set.remove(rec.name)
        orf_seq = Seq('')
        if len(read_hits) > 1:
            # If multiple orfs in one read then it's probably a frame
            # shift error or pseudo gene
            continue
        for index, hit in read_hits.sort('start').iterrows():
            # This for-loop is completely uneccesary thanks to the if statement
            # above.
            start = int(hit['start'])
            stop = int(hit['stop'])
            if start < stop:
                orf_seq += rec.seq[(start - 1):stop]
            elif start > stop:
                orf_seq = rec.seq[(stop - 1):start].reverse_complement() + orf_seq
            orf_rec = SeqRecord(id=hit['orf'], seq=orf_seq, description='')
            write(orf_rec, sys.stdout, 'fasta')
