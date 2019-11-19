import sys
import gffutils
#import pyfaidx
from Bio import SeqIO
from HTSeq import GFF_Reader
import textwrap
from ipdb import set_trace

genome = sys.argv[1]
gff = sys.argv[2]
cds_file = sys.argv[3]
prot_file = sys.argv[4]

# Clean exonerate output
gff_str=""
gene_ct = 1
for feat in GFF_Reader(gff):
    if feat.type == "gene":
        query_name = feat.attr["sequence"].strip()
        gene = f"{query_name}_{str(gene_ct)}"
        gene_ct += 1
        name = gene
        count = 1
    elif feat.type == "cds":
        name = f"{gene}.{count}"
        count += 1
    chrom = feat.iv.chrom
    start = feat.iv.start + 1
    end = feat.iv.end
    #phase = feat.get_gff_line().split("\t")[7]
    phase = feat.frame
    endings = dict(gene="\n",cds=f";Parent={gene}\n")
    if feat.type in ["gene", "cds"]:
        ending = endings[feat.type]
        gff_str += f"{chrom}\texonerate\t{feat.type}\t{start}\t{end}\t.\t{feat.iv.strand}\t{phase}\tID={name}{ending}"

parsed_gff_file = gff.replace(".gff", "_parsed.gff")
with open(parsed_gff_file, "w") as gff_o:
    gff_o.write(gff_str)

db = gffutils.create_db(parsed_gff_file, dbfn=":memory:")
fasta = genome
cds_seq = ""
cds_d = {}
gene = ""
for gene in db.features_of_type("gene", order_by="start"):
    gene_name = gene.id
    cds_d[gene_name] = ""
    for cds in db.children(gene_name, featuretype="cds"):
        cds_d[gene_name] += cds.sequence(genome)

# Save CDS sequences
with open(cds_file, "w") as cds_f:
    for gene, seq in cds_d.items():
        cds_f.write(f">{gene}\n")
        wrapped_seq = textwrap.fill(seq, 70) + "\n"
        cds_f.write(wrapped_seq)

# Translate and save protein sequences
prots = []
for seq in SeqIO.parse(cds_file, "fasta"):
    prot_s = seq.translate()
    prot_s.id = seq.id
    prot_s.description = ""
    prots.append(prot_s)

SeqIO.write(prots, prot_file, "fasta")
