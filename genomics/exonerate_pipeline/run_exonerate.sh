#!/usr/bin/env bash

# This script runs the analysis
# Usage:
#bash run.sh query_sequence.fa target_genome.fa

query=$1
genome=$2

# Run exonerate
exonerate --model protein2genome -q $query -t $genome --showalignment no --showvulgar no --showtargetgff yes > exonerate_out
