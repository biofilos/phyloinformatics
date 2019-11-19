import sys
import gffutils
#import pyfaidx
from Bio import SeqIO
from HTSeq import GFF_Reader

genome = sys.argv[1]
gff = sys.argv[2]
cds = sys.argv[3]
prot = sys.argv[4]

# Clean exonerate output
gff_str=""
gene_ct = 1
for feat in GFF_Reader("meh.gff"):
    if feat.type == "gene":
        gene = "gene" + str(gene_ct)
        gene_ct += 1
        name = gene
        count = 1
    elif feat.type == "cds":
        name = f"{gene}_{count}"
        count += 1
    chrom = feat.iv.chrom
    start = feat.iv.start + 1
    end = feat.iv.end
    phase = feat.get_gff_line().split("\t")[7]
    endings = dict(gene="\n",cds=f";Parent={gene}\n")
    if feat.type in ["gene", "cds"]:
        ending = endings[feat.type]
        gff_str += f"{chrom}\texonerate\t{feat.type}\t{start}\t{end}\t.\t{feat.iv.strand}\t{phase}\tID={name}{ending}"

db = gffutils.create_db(gff_str,dbfn=":memory:",from_string=True)
#fasta = pyfaidx.Fasta(genome)
fasta = genome
cds_seq = ""
cds_d = {}
gene = ""
for gene in db.features_of_type("gene", order_by="start"):
    gene_name = gene.id
    cds_d[gene_name] = ""
    for cds in db.children(gene_name, featuretype="cds"):
        cds_d[gene_name] += cds.sequence(genome)
