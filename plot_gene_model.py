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
    tempfile = "draft.svg"
    surface = cairo.SVGSurface(tempfile, 1280, 200)
    cr = cairo.Context(surface)
    cr.select_font_face(font, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    cr.set_font_size(fontsize)
    xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(text)
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
    features = ["CDS"]
    genes = {}
    feature_limits = dict(gff_type=features + ["gene", "transcript"])
    for rec in GFF.parse(gff, limit_info=feature_limits):
        for gene in rec.features:
            for transcript in gene.sub_features:
                transcript_start = transcript.location.start
                if transcript.sub_features:
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

                genes[transcript.id] = parsed_exons

    return genes


def parse_txt_css(css):
    font_params = {}
    parse = False
    for line in open(css):
        if parse:
            if "font-size" or "font-family" in line:
                param, value = line.strip().replace(";", "").split(": ")
                font_params[param] = value
            if len(font_params) ==2:
                break
        if ".name" in line:
            parse = True
    return font_params


def plot_gene(gene, coords, plot_obj, norm, gene_height, text_space, font_params):
    exon_height = gene_height
    fig_width = plot_obj.attribs["width"]
    middle_line = exon_height / 2
    gene_group = plot_obj.add(plot_obj.g(id="gene"))
    gene_len = max([x[-1] for x in coords]) * norm
    gene_group.add(plot_obj.line(start=(0,middle_line), end=(gene_len, middle_line), class_="line"))
    # Normalize the sizes so that they fit in the image
    last_end_coord = coords[-1][1]
    #norm = fig_width / last_end_coord
    for start, end in coords:
        gene_group.add(plot_obj.rect(insert=(start * norm, 0), class_="exon",
                                 size=(end * norm - start * norm, exon_height)))
    font = font_params["font-family"]
    font_size = int(font_params["font-size"].replace("px", ""))
    txt_w = get_text_metrics(font, font_size, gene)
    txt_y = middle_line #+ (txt_h / 2)
    txt = plot_obj.text("", insert=(end * norm + text_space, txt_y))
    txt.add(plot_obj.tspan(gene, dy=["0.35em"], class_="name"))
    gene_group.add(txt)


def plot_genes(genes_d, css_file, outfile, fig_width, gene_height, gene_space=10, text_space=2):
    total_height = gene_height * len(genes_d) + gene_space * (len(genes_d) - 1)
    font_params = parse_txt_css(css_file)
    font = font_params["font-family"]
    font_size = int(font_params["font-size"].replace("px", ""))
    max_txt_len = max([get_text_metrics(font, font_size, x) for x in genes_d.keys()])
    total_width = fig_width + max_txt_len + text_space
    plot = svgwrite.Drawing(outfile, size=(total_width, total_height))
    css = open(css_file).read()
    plot.embed_stylesheet(css)
    longest_gene = max([x[-1][1] for x in genes_d.values()])
    norm = fig_width / longest_gene
    for gene, coords in genes_d.items():
        plot_gene(gene, coords, plot, norm, gene_height, text_space, font_params)
    # Organize genes
    moved_space = 0
    for ix, gene in enumerate(plot.elements[1:]):
        gene.translate(0, moved_space)
        moved_space += (gene_height + gene_space)
    plot.save()

gff = argv[1]
width = int(argv[2])
height = int(argv[3])
svg_out = argv[4]

genes = extract_genes(gff)
plot_genes(genes, "style.css", svg_out, width, height)
print("DONE")
