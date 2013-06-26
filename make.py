#! /usr/bin/env python
"""Carry out a full analysis based on parameters defined in metadata.

"""

import os
from itertools import product
from pymake import Rule, maker

MIN_ORF_SIZE = 30
E_VALUE_CUTOFF = 1e-5

def read_col(path, col=0, sep='\t', skip=[]):
    """Return an iterable of the items in the specified column.

    """
    with open(path) as f:
        for i, line in enumerate(f):
            if i not in skip:
                yield line.split(sep)[col].strip()

samples = list(read_col('meta/samples.tsv', col=0, skip=[0]))
genes = list(read_col('meta/genes.list'))
refs = list(read_col('meta/refs.list'))
models = genes + refs
hit_files = ["hit/{sample}.orf{size}.{model}.hmmsearch_out".\
                 format(sample=s, size=MIN_ORF_SIZE, model=m)
             for s, m in product(samples, models)]
hit_files_str = " ".join(hit_files)
hit_seqs_files = ["seq/hit/fn/{sample}.orf{size}.{model}.fn".\
                 format(sample=s, size=MIN_ORF_SIZE, model=m)
             for s, m in product(samples, models)]


rules = [Rule("all",
              preqs=["res/hits_figures.png", "res/summary.tsv"]),
         Rule("res/hits_figures.png",
              preqs=["bin/make_hits_figures.py", "meta/samples.tsv",
                     "res/normtally.tsv"],
              recipe="{preqs[0]} -s {preqs[1]} -o {trgt} {preqs[2]}"),
         Rule("hit_seqs", preqs=hit_seqs_files),
         Rule(r"seq/hit/fn/(([^.]*)(\..*)?\.([^.]*))\.fn",
              preqs=["bin/fetch_frags.py", "hit/{0}.hmmsearch_out",
                     "seq/raw/{1}.fq", "seq/hit/fn/"],
              recipe=("{preqs[0]} -f fastq -e {E_VALUE_CUTOFF} "
                      "{preqs[1]} {preqs[2]} > {trgt}"),
              E_VALUE_CUTOFF=E_VALUE_CUTOFF),
         Rule(r"res/normtally\.tsv",
              preqs=["bin/data_compile.py", "hmm/components.tsv",
                     "hmm/models.tsv", "hmm/properties.tsv", "hmmsearch", "res/"],
              recipe=("{preqs[0]} -t normtally -e {E_VALUE_CUTOFF} "
                      "{preqs[1]} {preqs[2]} {preqs[3]} {hit_files} > {trgt}"),
              hit_files=hit_files_str, E_VALUE_CUTOFF=E_VALUE_CUTOFF),
         Rule(r"res/summary\.tsv",
              preqs=["bin/data_compile.py", "hmm/components.tsv",
                     "hmm/models.tsv", "hmm/properties.tsv",
                     "meta/samples.tsv", "hmmsearch", "res/"],
              recipe=("{preqs[0]} -t summary -s {preqs[4]} -e {E_VALUE_CUTOFF} "
                      "{preqs[1]} {preqs[2]} {preqs[3]} {hit_files} > {trgt}"),
              hit_files=hit_files_str, E_VALUE_CUTOFF=E_VALUE_CUTOFF),
         Rule("hmmsearch",
              preqs=hit_files),
         Rule(r"hit/([^.]*(\..*)?)\.([^.]*)\.hmmsearch_out",
              preqs=["bin/hmmsearch.sh", "hmm/{2}.hmm",
                     "seq/orf/{0}.fa", "hit/"],
              recipe=("{preqs[0]} {preqs[1]} {preqs[2]} {trgt}")),
         Rule(r"seq/orf/(([^.]*)(\..*)?)\.orf([0-9]*)\.fa",
              preqs=["bin/orfs.sh", "seq/raw/{0}.fq", "seq/orf/"],
              recipe=("{preqs[0]} -l {3} -f fastq < {preqs[1]} > {trgt}")),
         Rule(".*/",
              recipe="mkdir -p {trgt}",
              order_only=True),
         ]

if __name__ == "__main__":
    maker(rules)
