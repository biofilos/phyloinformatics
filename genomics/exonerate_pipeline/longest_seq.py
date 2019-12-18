from collections import defaultdict
from HTSeq import GFF_Reader
from sys import argv
from Bio import SeqIO

"""
Extracts the longest gene of every overlapping set of genes from a gff file
usage:
    python longest_seq.py gff_in cds_in prot_in gff_out cds_out prot_out
"""

def genes_overlap(gene1, gene2):
    """
    tests if two genes have overlapping coordinates
    :param gene1: HTSeq record
    :param gene2: HTSeq record
    :return: bool
    """
    overlap = False
    # Process genes if they are in the same chromosome
    if gene1.iv.chrom == gene2.iv.chrom:
        g1_coord = gene1.iv
        g2_coord = gene2.iv
        # This assumes that the GFF is sorted

        # If the start coordinate of the second gene is
        # smaller than the end coordinate of the first gene,
        # they overlap
        if g2_coord.start <= g1_coord.end:
            overlap = True
    return overlap


def get_longest(gene1, gene2):
    """
    Return the longest of two genes
    :param gene1: HTSeq record
    :param gene2: HTSeq record
    :return: HTSeq feature
    """
    if gene2.iv.length > gene1.iv.length:
        return gene2
    else:
        return gene1


# Set input/output files
gff_in = argv[1]
cds_in = argv[2]
aas_in = argv[3]
gff_out = argv[4]
cds_out = argv[5]
prot_out = argv[6]

# Initialize list of non-overlapping genes
gene_d = defaultdict(list)
gene_list = []
all_genes = []
for rec in GFF_Reader(gff_in):
    if rec.type == "gene":
        all_genes.append(rec)
        gene_d[rec.iv.chrom].append(rec)

# Sort genes in each chromosome
for chrom in gene_d:
    gene_d[chrom] = sorted(gene_d[chrom],key=lambda x: x.iv.start)
for genes in gene_d.values():
    for rec in genes:
        if not gene_list:
            gene_list.append(rec)
        else:
            prev = gene_list[-1]
            curr = rec
            if genes_overlap(prev, curr):
                gene_list[-1] = get_longest(prev, curr)
            else:
                gene_list.append(curr)

# Set sequence dictionaries
cds_d = SeqIO.to_dict(SeqIO.parse(cds_in, "fasta"))
prot_d = SeqIO.to_dict(SeqIO.parse(aas_in, "fasta"))
# Write new GFF
with open(gff_out, "w") as gff_o, open(cds_out, "w") as cds_o, open(prot_out, "w") as prot_o:
    for rec in gene_list:
        gff_line = rec.get_gff_line()
        gff_line = gff_line.replace(' "', "=").replace('"', "")
        gff_o.write(gff_line)
        cds_seq = cds_d[rec.name]
        prot_seq = prot_d[rec.name]
        SeqIO.write(cds_seq, cds_o, "fasta")
        SeqIO.write(prot_seq, prot_o, "fasta")
# Write new sequence files
