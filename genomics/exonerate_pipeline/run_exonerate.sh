#!/usr/bin/env bash

# This script runs the analysis
# Usage:
#bash run_exonerate.sh query_sequence.fa target_genome.fa cds_file protein_file

query=$1
genome=$2
cds_file=$3
prot_file=$4
exonerate_out=exonerate.txt
# Run exonerate
exonerate --model protein2genome -q $query -t $genome --showalignment no --showvulgar no --showtargetgff yes > $exonerate_out
exonerate_gff=`basename $exonerate_out .txt`.gff
#cat $exonerate_out | grep -v Command | grep -v Hostname  | grep -P -v "\tsimilarity\t" | grep -v "#" | grep -v splice > $exonerate_gff
cat $exonerate_out | grep -v Command | grep -v Hostname  | grep -v "\-\- completed exonerate analysis" > $exonerate_gff

# Parse exonerate output and generate sequence files
python extract_seqs.py $genome $exonerate_gff $cds_file $prot_file
