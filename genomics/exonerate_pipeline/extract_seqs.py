import sys
import gffutils
from Bio import SeqIO
from HTSeq import GFF_Reader
import textwrap


"""
The purpose of this script is to further format the pre-processed GFF file obtained
from the run_exonerate.sh script, and extract cds and protein sequences from its
results
Usage:
    python extract_seqs.py genome.fa pre-processed_exonerate_out.gff out_cds.fa 
    out_proteins.fa out_gff.gff
"""

# 1. Load files
genome = sys.argv[1]
gff = sys.argv[2]
cds_file = sys.argv[3]
prot_file = sys.argv[4]
final_gff = sys.argv[5]

# 2. Clean exonerate output
# Initialize gff text
gff_str=""
# Initialize gene counter
gene_ct = 1
# Read each record (assume the file is ordered, so the cds features below each
# gene feature will be interpreted as belonging to such gene feature
for feat in GFF_Reader(gff):
    if feat.type == "gene":
        # If the feature is a gene, extract and format its name

        # This line extracts the id of the query sequence that generated this
        # gene annotation
        query_name = feat.attr["sequence"].strip()
        # Set the query id as the gene, plus an underscore and a consequtive
        # number (to make all the gene names unique)
        gene = f"{query_name}_{str(gene_ct)}"
        gene_ct += 1
        name = gene
        count = 1
    elif feat.type == "cds":
        # Set the name of the cds in each gene as the name of the gene,
        # followed by a dot ".", followed by a consequtive number
        name = f"{gene}.{count}"
        count += 1
    # Extract coordinates of the feature
    chrom = feat.iv.chrom
    # Add one to the start coordinate, because in this library, the start
    # coordinate is zero-based, and in GFF, they are 1-based
    start = feat.iv.start + 1
    end = feat.iv.end
    phase = feat.frame
    # Set attribute options for genes and cds (see below)
    endings = dict(gene="\n",cds=f";Parent={gene}\n")
    if feat.type in ["gene", "cds"]:
        # If the feature is a gene, do not include any attributes (other than
        # the ID). If the feature is a cds, add the "Parent" attribute, to link
        # said cds to its gene
        ending = endings[feat.type]
        # Assemble gff line
        gff_str += f"{chrom}\texonerate\t{feat.type}\t{start}\t{end}\t.\t{feat.iv.strand}\t{phase}\tID={name}{ending}"

# 3. Write formatted GFF file
with open(final_gff, "w") as gff_o:
    gff_o.write(gff_str)

# 4. Load parsed GFF as a gffutils database
db = gffutils.create_db(final_gff, dbfn=":memory:")
# Set variables
fasta = genome
cds_seq = ""
cds_d = {}
gene = ""
# 5. Go through the genes in the GFF
for gene in db.features_of_type("gene", order_by="start"):
    # Obtain gene name
    gene_name = gene.id
    # Initialize data for the gene with empty text
    cds_d[gene_name] = ""
    for cds in db.children(gene_name, featuretype="cds"):
        # Add the dna sequence of each cds to the text of each gene
        cds_d[gene_name] += cds.sequence(genome)

# 6. Save CDS sequences
with open(cds_file, "w") as cds_f:
    for gene, seq in cds_d.items():
        # Set the name of the gene (from the GFF) as the id for the resulting
        # fasta files
        cds_f.write(f">{gene}\n")
        # Format sequences text
        wrapped_seq = textwrap.fill(seq, 70) + "\n"
        # Write sequence to file
        cds_f.write(wrapped_seq)

# 7. Translate and save protein sequences
prots = []
for seq in SeqIO.parse(cds_file, "fasta"):
    # Translate sequence
    prot_s = seq.translate()
    # Format id of each sequence, so that thay match the cds file (useful if
    # performing dN/dS
    prot_s.id = seq.id
    prot_s.description = ""
    prots.append(prot_s)
# Write protein sequences
SeqIO.write(prots, prot_file, "fasta")
