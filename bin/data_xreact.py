#!/usr/bin/env python
"""Determines if any models share hits.

This is usually a bad thing, except in the case where models are
meant to cross-react.

"""

import numpy as np
from itertools import permutations
from glob import glob
from data_compile import compile, forwhich
from data_compile import get_all_hmmsearch_out_paths
from data_compile import load_and_merge


def uniq_orf_ids(df):
    return df.seq_set + df.orf

hits, tally, norm_tally, results = compile(get_all_hmmsearch_out_paths('hit/'),
                                           1e-5,
                                           load_and_merge(glob("hmm/*.tsv")))
hits['uniq_orf_id'] = hits.sample + hits.orf

component_list = list(set(hits.component))
orfs = {}
for component in component_list:
    orfs[component] = set(forwhich(hits, component=component).uniq_orf_id)
matrix = np.empty(shape=(len(component_list), len(component_list)))
for i, component_i in enumerate(component_list):
    for j, component_j in enumerate(component_list[:i]):
        intersect_size = float(len(orfs[component_i] & orfs[component_j]))
        matrix[i, j] = intersect_size / len(orfs[component_i])
        matrix[j, i] = intersect_size / len(orfs[component_j])
xreactions = [(component_list[i], component_list[j], matrix[i, j])
              for i, j in permutations(range(len(component_list)), 2)
              if matrix[i, j] > 0]
for item in xreactions:
    print(item)
