#!/usr/bin/env bash

# This script runs the analysis
# Usage:

#bash run_exonerate.sh query_sequence.fa target_genome.fa exonerate_out cds_file protein_file gff_output

# This specifies the file names in the same order as they were
# typed in the command (see usage)
query=$1
genome=$2
exonerate_out=$3
cds_file=$4
prot_file=$5
gff_out=$6

root_path=$(cd `dirname $0` && pwd)
# Run exonerate

# Exonerate extra options
# Modify the following variable with more options from exonerate. NOTE: Only change
# what is inside the double quotes
exonerate_options=""

# Exonerate comes with a variety of models to be run. Here, the model is self
# explanatory. More information can be found at https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide
exonerate --model protein2genome $exonerate_options -q $query -t $genome --showalignment no --showvulgar no --showtargetgff yes > $exonerate_out
# Replace the extension of the raw exonerate output, and append .gff. This
# will be the name of the final GFF
exonerate_gff=$(echo "$exonerate_out" | cut -f 1 -d '.').gff

# Parse relevant lines from GFF information, and save them 
cat $exonerate_out | grep -v Command | grep -v Hostname  | grep -v "\-\- completed exonerate analysis" > $exonerate_gff
# Parse exonerate output and generate sequence files
python ${root_path}/extract_seqs.py $genome $exonerate_gff $cds_file $prot_file $gff_out 
