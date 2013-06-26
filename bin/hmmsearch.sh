#! /usr/bin/env bash
# Run a hmmsearch and mung the results.

# TODO: --tblout should be --pfamout now
hmmsearch --tblout $3 --tformat FASTA --cpu 2 $1 $2

sed -i '/^#/d' $3
