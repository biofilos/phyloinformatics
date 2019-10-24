# phyloinformatics
Genomics, phylogenetics, and data visualization code developed at the Laboratory for Phyloinformatics at RIKEN BDR

## visualization
### plot_gene_model.py
Description  
Plots exon (using CDS features) structure of the genes in a GFF file and saves a svg file. The GFF file should contain the features "gene", "transcript", and "CDS" (see sample.gff). The text next to each gene structure will be the transcript ID from the GFF file.  
Dependencies: BCBio, svgwrite, cairo (and pycairo)  
Usage: `pyton plot_gene_model.py gff_file width_pixels height_pixels svg_file.svg`  
Styling: A CSS file named `style.css` is included. In it, basic styling parameters can be adjusted for the text, exon, and line that spans each gene with the classes ".name", ".exon", and ".line", respectively. With it, relevant colors and width of elements can be customized.
