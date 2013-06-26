#! /usr/bin/env bash
# Takes a single integer argument, the minimum orf length.
# FASTQ formatted sequences are taken to stdin and FASTA formatted ORFs are
# output to stdout.

# getorf requires libimf.so
# Or you'll get:
#
#   getorf: error while loading shared libraries: libimf.so: cannot open shared
#   object file: No such file or directory

# On flux this is accessed with `module load intel-comp`.

# Many FASTA formatted files will need to be considered "pearson" format instead in order
# to not mangle the ID.  Use the '-f' option.
# See http://emboss.sourceforge.net/docs/themes/SequenceFormats.html

MIN_ORF_LEN=30

while getopts ":l:f:" opt; do
	case $opt in
		l) MIN_ORF_LEN=$OPTARG ;;
		f) FORMAT=$OPTARG ;;
		*) echo "Invalid option: -$OPTARG" >&2
		   exit 1 ;;
	esac
done

if [ ! -z "$FORMAT" ]; then
	FORMAT_OPTION="-sformat1 $FORMAT"
fi

getorf -filter $FORMAT_OPTION -minsize $MIN_ORF_LEN < /dev/stdin | \
sed -r '/^>/s:_[0-9]+ \[([0-9]+) - ([0-9]+)\].*:(\1-\2):g' > /dev/stdout
