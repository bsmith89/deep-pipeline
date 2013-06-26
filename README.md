An example bioinformatics pipeline using pymake.

## Dependencies

  +  PyMake (github.com/bsmith89/pymake)
  +  HMMER3 (http://hmmer.janelia.org/software)
  +  getorf (an EMBOSS tool: http://emboss.sourceforge.net/download/#Stable/)
  +  python3 (and preferrably a virtualenv)
       +  Pandas
       +  Biopython

## Quickstart

Once you have all of those tools installed,
move sequence files into 'seq/raw/', hmm's and
all of the hmm metadata tsv's into 'hmm/' and
edit the sample metadata, and the reference
and target gene lists in 'meta/' if so desired.

Then:

	./make.py -n res/summary.tsv

in the root directory to get a printout of the steps that will be carried
to create a summary of the target/reference gene hits in the samples.

This command should not give you any errors, although the lack of errors does
*not* mean that the pipeline will complete successfully.

When you are satisfied that all of the necessary utilities are installed,
and that the commands printed with the '-n' flag set (meaning: 'dry run')
will not be destructive (the pipeline *will* overwrite extant files, although
this may be changed in future versions of PyMake) do:

	./make.py res/summary.tsv


## Files have a consistant naming scheme:

Any file that is downstream of a sequence file
(whether it is a sequence file itself or not)
will be named

	<SAMPLE>[.<FILTER1>[.<FILTER2>[...]]].<EXTENSION>

for instance

	SAMPLE1.fn

would be the raw, unmodified nucleotide sequnce (unaligned) of
SAMPLE1's reads.

`<FILTERS>` can be a variety of things, for instance

    SAMPLE1.fused.fq

would be the raw, paired-end reads from SAMPLE1 fused together
across their overlap to make a FASTQ formatted file (nucleotide).

	'SAMPLE2.fused.orf30.fa'

would be the fused reads from SAMPLE2, run through an ORF finder with
a minimum length of 30 resulting in a FASTA formatted file of amino-acid
sequences. (The ORF sequences, as it turns out, will have the same read
names as their source reads, but with '(START-END)' appended to the end to
indicate the coordinates in nucleotide space which defines the ORF.

	'SAMPLE1.fused.orf30.TIGR02376.hmmsearch_out'

is the hmmsearch tabular output file (with all '# Comments' removed) for
SAMPLE1, processed as described above and searched against the TIGR02376
profile-HMM.

	'SAMPLE1.fused.orf30.TIGR02376.fn'

Are the nucleotide sequences of the ORFs which hit TIGR02376 from SAMPLE1.

As you can see, we can parse out a "middle" portion which is just the
filtering, and an optional penultimate item, which is the model, followed
by the extension which indicates what type of data is contained in the
file.

By convention and coding, `<SAMPLE>` should not contain any periods, and nor should
`<MODEL>`.
