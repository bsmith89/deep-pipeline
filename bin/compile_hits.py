#!/usr/bin/env python
"""Extract and compile hmmsearch results with metadata.

First step in analyzing hmmsearch data.

"""

import os
import os.path
import sys
import optparse
import pandas as pd
import numpy as np
from itertools import product
from glob import glob
from warnings import warn

_DEFAULT_E_VALUE_CUTOFF = 1e-5


def forwhich(df, **kwargs):
    """Return the rows of DataFrame df which satisfy the conditions.

    Conditions are passed as keyword arguments to the function.  Only
    equality conditions work.
    """
    return df[np.all(np.row_stack([df.get(key) == value
                                   for key, value in kwargs.items()]), 0)]


def compile(hit_file_paths, e_value_cutoff, model_details, samples_data=None):
    """Compile the hmmsearch results in *hit_file_paths*.

    If a DataFrame object is passed as *meta*,
    it is merged against the hits table.

    """
    samples = set()
    models = set()
    crosses = set()
    hits_tables = []
    for path in hit_file_paths:
        details = os.path.basename(path).split('.')
        sample, *processing, model, extension = details
        samples.add(sample)
        models.add(model)
        crosses.add((sample, model))
        if os.path.getsize(path) == 0:
            continue
        print("Opening %s" % path, file=sys.stderr)
        table = pd.read_table(path, sep='\s*', header=None)
        hits = table.select(lambda col: col in [0, 2, 4, 5, 6], axis=1)
        hits.columns = ['orf', 'model', 'e_value', 'score', 'bias']
        hits['sample'] = sample
        hits['processing'] = ','.join(processing)
        hits_tables.append(hits)
    missing_crosses = set(product(samples, models)) - crosses
    if missing_crosses:
        warn(UserWarning(("Incomplete Data: Not all samples-by-models "
                          "combinations were provided.  Specifically, we're missing:\n"
                          "{}").format(missing_crosses)))
    hits = pd.concat(hits_tables, ignore_index=True)
    hits = hits[hits.e_value <= e_value_cutoff]
    hits['read'], hits['start'], hits['stop'] = \
        np.array(hits.orf.str.\
            match('^(.*)\(([0-9]*)-([0-9]*)\)$').values.tolist()).T
    hits = pd.merge(hits, model_details[['model', 'component']].drop_duplicates())
    # Tally uniq hits to all of the models which point at a given component.
    tally = hits.groupby(['component', 'sample']).\
        apply(lambda df: len(set(df.sample + df.orf.str.replace('_[0-9]+', '')))).unstack()
    tally = tally.fillna(value=0)
    norm_tally = tally / tally[tally.index.isin(
        model_details[model_details.property ==
                      'bjs_bact_ubiq_1cpy'].component)].mean()
    results = pd.concat(dict(tally=tally.unstack(),
                            norm_tally=norm_tally.unstack()),
                        axis=1).reset_index()
    out = [hits, tally.reset_index(), norm_tally.reset_index(), results]
    if samples_data:
        treat = pd.merge(results, samples_data)\
            [['treat', 'sample', 'component', 'norm_tally', 'tally']].groupby(('treat', 'component'))
        treat_mean = treat.mean()
        treat_std = treat.std()
        summary = pd.merge(treat_mean, treat_std,
                           left_index=True, right_index=True,
                           suffixes=['_mean', '_std']).reset_index()
        out.append(summary)
    return tuple(out)


def main():
    table_choices = ["tally", "normtally", "both", "allhits", "summary"]
    DEFAULT_TABLE_OUT = "summary"
    usage = ("usage: %prog [options] COMPONENTS.tsv MODELS.tsv PROPERTIES.tsv "
             "HIT_FILE1 [HIT_FILE2 [HIT_FILE3 ...]]")
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-e", "--e-value-cutoff", dest="e_value_cutoff",
                      default=_DEFAULT_E_VALUE_CUTOFF,
                      type='float',
                      help=("the maximum E-value considered a hit. "
                            "DEFAULT: {}").format(_DEFAULT_E_VALUE_CUTOFF))
    parser.add_option("-t", "--table-output", dest="table", type="choice",
                      choices=table_choices, default=DEFAULT_TABLE_OUT,
                      help=("which results table to write. "
                            "Choose from {0}. "
                            "DEFAULT: '{1}'").format(table_choices, DEFAULT_TABLE_OUT))
    parser.add_option("-s", "--sample-sheet", dest="samples_path",
                      help=("path to the sample metadata. "
                            "Only required for '--table=summary' "))
    opts, args = parser.parse_args()

    components_path, models_path, properties_path, *hit_file_paths = args
    model_details = pd.merge(pd.merge(pd.read_table(components_path),
                                      pd.read_table(models_path)),
                             pd.read_table(properties_path))
    compile_args = [hit_file_paths, opts.e_value_cutoff, model_details]
    if opts.table == 'summary':
        compile_args += [pd.read_table(opts.samples_path)]
    hits, tally, norm_tally, results, *summary = compile(*compile_args)
    if summary:
        summary = summary[0]  # Since '*summary' unpacked value collection yields a list.
    results_table = dict(zip(table_choices, [tally, norm_tally, results, hits, summary]))
    results_table[opts.table].to_csv(sys.stdout, index=False, na_rep='NA', sep='\t')

if __name__ == '__main__':
    main()
