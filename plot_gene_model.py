"""
Given a GFF file, plot the genes in it in a SVG file
Optional: specify colors for gene features (exons, introns, etc)
in a JSON file
Usage:
python plot_gene_model.py <gff_file> <json_color_specs optional>

Notes:
    The GFF must have at least "gene" and "transcript" types (column 3)
"""

from BCBio import GFF
from sys import argv
import json
import os
from ipdb import set_trace
import svgwrite
import tkinter as Tkinter

try:
    import tkFont
except ImportError:
    import tkinter.font as tkFont

def get_text_metrics(family, size, text):
    # initialize Tk so that font metrics will work
    tk_root = Tkinter.Tk()
    font = None
    key = (family, size)
    font = tkFont.Font(family=family, size=size)
    assert font is not None
    (w, h) = (font.measure(text), font.metrics('linespace'))
    return (w, h)
def extract_genes(gff, colors_d):
    """
    Extract gene models from gff files. If a color_specs file
    is provided, include the features specified in the
    color_specs file
    :param gff: path to gff file
    :param colors_d: dictionary of color specifications: {"feature": "color"}
    :return: list of BCBio genes
    """
    if colors_d:
        features = color_specs.keys()
    else:
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


def svg(in_coords, out_file, fig_width, fig_height):
    plot = svgwrite.Drawing(out_file, size=(fig_width, fig_height))
    css = "{font-size:14px;font-family:'sans-serif'}"
    exon_height = fig_height
    gene_group = plot.add(plot.g(id="gene"))
    # Normalize the sizes so that they fit in the image
    last_end_coord = in_coords[-1][1]
    norm = fig_width / last_end_coord
    middle_line = exon_height / 2
    #line_style = "stroke:#000000;stroke-opacity:1"
    plot.add(plot.line(start=(0,middle_line), end=(fig_width, middle_line),stroke_width=50,stroke_opacity=1,stroke="black"))
    for start, end in in_coords:
        gene_group.add(plot.rect(insert=(start * norm, 0),
                                 size=(end * norm - start * norm, exon_height),fill="coral"))
    plot.save(pretty=True)


gff = argv[1]

if len(argv) == 3 and os.path.exists(argv[2]):
    json_path = argv[2]
    color_specs = json.load(open(json_path))
else:
    color_specs = False

genes = extract_genes(gff, color_specs)
svg(genes["Chipu0000003.t1"], "delete.svg", 1000, 300)
print("DONE")
