"""
Given a GFF file, plot the genes in it in a SVG file
Optional: specify colors for gene features (exons, introns, etc)
in a JSON file
Usage:
python plot_gene_model.py <gff_file> width height <out_file.svg>

Notes:
    The GFF must have at least "gene" and "transcript" types (column 3)
"""

from BCBio import GFF
from sys import argv
#import json
import os
from ipdb import set_trace
import svgwrite
#import tkinter as Tkinter
import cairo
#try:
#    import tkFont
#except ImportError:
#    import tkinter.font as tkFont

##def get_text_metrics(family, size, text):
##    # initialize Tk so that font metrics will work
##    tk_root = Tkinter.Tk()
##    font = None
##    key = (family, size)
##    font = tkFont.Font(family=family, size=size)
##    assert font is not None
##    (w, h) = (font.measure(text), font.metrics('linespace'))
##    return w
def get_text_metrics(font, fontsize, text):
    """
    Measures the width in pixels of a strting of text
    :param font: Font name
    :param fontsize: Font size (int)
    :param: text: Text string to be measured
    """
    # This function was taken from https://stackoverflow.com/questions/24337531/how-to-determine-text-width-and-height-when-using-svgwrite-for-python/48387748
    # Cairo needs to write a file to be measured (this will be removed)
    tempfile = "draft.svg"
    # Generate surface
    surface = cairo.SVGSurface(tempfile, 1280, 200)
    cr = cairo.Context(surface)
    cr.select_font_face(font, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    cr.set_font_size(fontsize)
    xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(text)
    # Remove temporary file
    os.remove(tempfile)
    return width

def extract_genes(gff):
    """
    Extract gene models from gff files. If a color_specs file
    is provided, include the features specified in the
    color_specs file
    :param gff: path to gff file
    :return: list of BCBio genes
    """
    # Feature to be treated as exons
    features = ["CDS"]
    # Initialize genes dictionary 
    genes = {}
    # Only parse features in this list (plus CDS)
    feature_limits = dict(gff_type=features + ["gene", "transcript"])
    # Parse each record in the GFF
    for rec in GFF.parse(gff, limit_info=feature_limits):
        # Go through each gene in the record
        for gene in rec.features:
            # Process each transcript
            for transcript in gene.sub_features:
                transcript_start = transcript.location.start
                # Parse exons
                if transcript.sub_features:
                    # Capture exon coordinates. Subtract the start coordinates
                    # of the transcript, so that they start with zero
                    exons = [(x.location.start.real - transcript_start,
                              x.location.end.real - transcript_start) for x in transcript.sub_features
                             if x.type in features]
                else:
                    # If there are no exons, consider the start and end of the transcript as the coordinates
                    # of the only exon
                    exons = [(transcript.location.start - transcript_start,
                              transcript.location.end - transcript_start)]
                # Here I will have to put code to parse out UTRs
                # for now, make the exons start with 0
                sorted_exons = sorted(exons, key=lambda x: x[0])
                first_start = sorted_exons[0][0]
                if first_start > 0:
                    parsed_exons = [(x[0] - first_start,
                                     x[1] - first_start) for x in sorted_exons]
                else:
                    parsed_exons = sorted_exons
                # Add data to dictionay
                genes[transcript.id] = parsed_exons

    return genes


def parse_txt_css(css):
    """
    Extract font name and font size from style css
    :param css: Path to css file
    return dict()
    """
    # Initialize dictionary
    font_params = {}
    # start with parsing flag down
    parse = False
    # Read css file
    for line in open(css):
        # If parse flag is up (true), start parsing
        if parse:
            if "font-size" or "font-family" in line:
                param, value = line.strip().replace(";", "").split(": ")
                font_params[param] = value
            # If two values have been parsed (font-size and font-family), stop
            if len(font_params) ==2:
                break
        # Up the parsing flag when the file hits the .name class
        # which contains text information
        if ".name" in line:
            parse = True
    return font_params


def plot_gene(gene, coords, plot_obj, norm, gene_height, text_space, font_params):
    """
    Plots one gene in a svgwrite drawing
    :param gene: gene name
    :param coords: List of tupes with start and end positions of each exon
    :param plot_obj: svgwrite drawing object
    :param norm: normalization factor. This will be used to translate from genomic coordinates to plot coordinates (usually calculated by max(gene_lengths) / figure_width)
    :param gene_height: height (in pixels) of exons
    :param text_space: Space between the end of the gene and the start of the text
    :param font_params: dictionary containing font-family and font-size from css
    return None
    """
    # Extract figure width from drawing
    fig_width = plot_obj.attribs["width"]
    # Plot a horizontal line through the middle of the plot
    middle_line = gene_height / 2
    # Encapsulate shapes belonging to each gene in a group
    gene_group = plot_obj.add(plot_obj.g(id="gene"))
    # Get gene length (end coordinate of last exon)
    gene_len = max([x[-1] for x in coords]) * norm
    # Add line to gene group
    gene_group.add(plot_obj.line(start=(0,middle_line), end=(gene_len, middle_line), class_="line"))
    last_end_coord = coords[-1][1]
    # Process each exon
    for start, end in coords:
        # Adds exon to gene group
        gene_group.add(plot_obj.rect(insert=(start * norm, 0), class_="exon",
                                 size=(end * norm - start * norm, gene_height)))
    # Get font parameters and measure width of text
    font = font_params["font-family"]
    font_size = int(font_params["font-size"].replace("px", ""))
    txt_w = get_text_metrics(font, font_size, gene)
    # Vertical position of text at the middle line
    txt_y = middle_line
    # Create text object
    txt = plot_obj.text("", insert=(end * norm + text_space, txt_y))
    # Put text in Tsapn element
    # The dy parameter is included to approximate a text alignment in the
    # middle of the figure
    txt.add(plot_obj.tspan(gene, dy=["0.35em"], class_="name"))
    # Add text to gene group
    gene_group.add(txt)


def plot_genes(genes_d, css_file, outfile, fig_width, gene_height, gene_space=10, text_space=2):
    """
    Plot a set of genes and save resulting svg file
    :param genes_d: dictionary with gene coordinates (e.g.: {"gene1":[(start, end)...})
    :param css_file: css file to use for styling
    :param outfile: output svg file
    :param fig_width: figure width (pixels)
    :param gene_height: height of each gene (pixels)
    :param gene_space: distance between genes (pixels)
    :param text_space: distance between the end of each gene and its name (pixels)
    """
    # Calculate height of all genes
    all_genes_height = gene_height * len(genes_d)
    # Calculate height of all the spaces between genes
    all_spaces_height = gene_space * (len(genes_d) - 1)
    # Calculate figure length
    total_height = all_genes_height + all_spaces_height
    # Extract font parameters from css file
    font_params = parse_txt_css(css_file)
    font = font_params["font-family"]
    font_size = int(font_params["font-size"].replace("px", ""))
    # Calculate the width of the longest text (longest gene name)
    max_txt_len = max([get_text_metrics(font, font_size, x) for x in genes_d.keys()])
    # Calculate figure width
    total_width = fig_width + max_txt_len + text_space
    # Initialize drawing obejct
    plot = svgwrite.Drawing(outfile, size=(total_width, total_height))
    # Open and embed css in drawing
    css = open(css_file).read()
    plot.embed_stylesheet(css)
    # Calculate normalization factor to translate genomic coordinates to figure
    # coordinates
    longest_gene = max([x[-1][1] for x in genes_d.values()])
    norm = fig_width / longest_gene
    # Plot each gene
    for gene, coords in genes_d.items():
        plot_gene(gene, coords, plot, norm, gene_height, text_space, font_params)
    # Organize genes. At this poing, all genes have been plotted at the
    # beginning of the plot
    moved_space = 0
    # Move genes, so that they do not overlap
    for ix, gene in enumerate(plot.elements[1:]):
        gene.translate(0, moved_space)
        moved_space += (gene_height + gene_space)
    # Save image
    plot.save()

# Parse arguments
gff = argv[1]
width = int(argv[2])
height = int(argv[3])
svg_out = argv[4]
style = argv[0].replace("plot_gene_model.py","style.css")
# Parse GFF as dictionary of genomic coordinates
genes = extract_genes(gff)
# Plot genes
plot_genes(genes, style, svg_out, width, height)
# Finish
print(f"Figure saved: {svg_out}")
