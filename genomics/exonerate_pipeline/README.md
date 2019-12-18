## Exonerate pipeline

This short pipeline is useful to obtain potential homologous genes of a query protein(s) equence(s), in an unnanotated genome. For example, if one wants to build phylogenetic trees of a gene family, but in one of the species of interest, a genome sequence file (fasta) is the only file available (no proteome, cds, or GFF files).  

A conda environment is recommended to install all the required programs and libraries. A useful guide can be found [here](https://bioconda.github.io/user/install.html).  
### Requirements
* exonerate (tested with v2.4.0)
* Python 3
* HTSeq (Python library)
* gffutils (Python library)
* BioPython (Python library)

After creating a conda environment (plus the bioconda channels), all the requirements can be installed with: `conda install gffutils biopython htseq exonerate`.  

### Usage
NOTE: Remember to load the conda environment where the requirements were installed.  
The running script `run_exonerate.sh` will run exonerate and the python code required to clean the GFF and extract mRNA and protein sequences. An example run (input files included) follows:  
`bash run_exonerate.sh query.fa Mus_musculus.GRCm38.dna.chromosome.13.fa exonerate_out.txt cds.fa proteins.fa gff_out.gff`.  

#### Input:
* `query.fa`: Query **protein** sequence(s) you want to find potential homologous genes for.
* `genome.fa`: Fasta file of the **genomic sequences** where you want to find potential homologous genes in.
#### Output:
* `exonerate_out.txt`: Raw output from exonerate (only gff).
* `exonerate_out.gff`: Processed output from exonerate, including only relevant GFF data. This file is not specified by the user. The name will be the same as in the exonerate output (e.g. `exonerate_out.txt`, replacing its extension with .gff.
* `cds.fa`: **DNA** coding sequences of the potential homologous genes obtained in the target genome.
* `proteins.fa`: Translated **protein** sequences from `cds.fa`
* `gff_out.gff`: Parsed **GFF** file of the potential homologous genes obtained in the target genome.  
The ID attribute of the resulting GFF file follow the format:
* In gene features: \<sequence id from `query.fa`\>_\<number\>.
* In cds features: \<gene ID>.\<number\>

### Filtering out overlapping genes
Althought not common, exonerate can output overlapping genes. In those cases, those overlapping genes could be isoforms (though not necessarily), and the longest of them might be wanted. The script `longest_seq.py` included in this directory provides this functionality. `longest_seq.py` extracts the longest gene from each set of overlapping genes in a GFF, and writes GFF, CDS, and protein files with the non-overlapping genes. Genes that did not overlap are returned as well.  
Usage: `python longest_seq.py gff_in.gff cds_in.fa protein_in.fa filtered_gff.gff filtered_cds.fa filtered_proteins.fa`

### IMPORTANT
By design, the output from exonerate was not processed, only formatted. That means that the resulting coding and protein sequences should be further tested to make sure that they are likely homologs of the query sequence. As a suggestion, the gene legnths of both query and results can be inspected for big discrepancies. Also, domain annotation with Pfam can be performed for both query (and known homologs) and the results of this pipeline. Finally, a shallow phylogenetic analysis of query sequence (and known homologs) can be performed including the results from this pipeline, to test if the phylogenetic relationships between all the sequences are "reasonable".  
All that said, this is Biology, and well, Biology is complicated. There is no automation bioinformatics tool that can do the work of rigurous scientific reasoning. As such, the last step on this pipeline is your expertise and good judgement.

