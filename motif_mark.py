#!/usr/bin/env python

# Importing modules
# in bash terminal run "pip install pycairo"
import cairo
import re
import argparse

# Setting up argparse inputs
parser = argparse.ArgumentParser(description="'-fa' to specify fasta file; '-m' to specify motif mark file")
parser.add_argument("-fa", "--fasta_file", help="input file with fasta sequences", required = True)
parser.add_argument("-m", "--motif_mark_file", help="input file with motif marks", required = True)
args = parser.parse_args()

fasta_input = args.fasta_file
motif_input = args.motif_mark_file

###############################################################
################# BUILDING MOTIF CONVERSION FUNCTION BELOW

def motif_conversion(amotif):

    reg_building = ''

    for x in range(len(amotif)):
        if (len(vIUPAC[amotif[x]])) == 1:
            reg_building += (vIUPAC[amotif[x]][0])
        
        elif (len(vIUPAC[amotif[x]])) > 1:
            reg_building += "["
            for y in range(len(vIUPAC[amotif[x]])):
                reg_building += (vIUPAC[amotif[x]][y])
            reg_building += "]"

    return reg_building

########### END OF MOTIF CONVERSION FUNCTION
###############################################################
######## Start of exon_placer function

def exon_placer(aread):

    exon_regex = "[A-Z]+"

    test_exon_indices = [match.start() for match in re.finditer(exon_regex, aread)]
    test_exon = re.findall(exon_regex, aread)

    anexon_reg_found = re.search(exon_regex, aread)

    anexon = anexon_reg_found[0]

    exon_start_indices = [match.start() for match in re.finditer(exon_regex, aread)]
    exon_start = exon_start_indices[0]
    exon_end = exon_start + len(anexon)

    # list to output: the exon sequence, start position and end position 
    read_exon_start_end = [anexon, exon_start, exon_end]
    
    return read_exon_start_end


# so exon_placer fn() outputs a list of 3 elements:
# the exon found, the start position and the end position

###### End of exon_placer function
#################################################################


# Dict to translate ambiguous IUPAC nt symbols
# Credit to wikipedia: https://en.wikipedia.org/wiki/Nucleic_acid_notation)
# Although the dictionary was expanded to include lowercase nt symbols


vIUPAC = {
    "A":["A"            ],
    "C":[    "C"        ],
    "G":[        "G"    ],
    "T":[            "T"],
    "U":[            "U"],
    "W":["A",        "T"],
    "S":[    "C","G"    ],
    "M":["A","C"        ],
    "K":[        "G","T"],
    "R":["A",    "G",   ],
    "Y":[    "C",    "T"],
    "B":[    "C","G","T"],
    "D":["A",    "G","T"],
    "H":["A","C",    "T"],
    "V":["A","C","G",   ],
    "N":["A","C","G","T"],
    "Z":[               ],

    "a":["a"            ],
    "c":[    "c"        ],
    "g":[        "g"    ],
    "t":[            "t"],
    "u":[            "u"],
    "w":["a",        "t"],
    "s":[    "c","g"    ],
    "m":["a","c"        ],
    "k":[        "g","t"],
    "r":["a",    "g",   ],
    "y":[    "c",    "t"],
    "b":[    "c","g","t"],
    "d":["a",    "g","t"],
    "h":["a","c",    "t"],
    "v":["a","c","g",   ],
    "n":["a","c","g","t"],
    "z":[               ],
}

########################################################################
# Building list of motifs from input file
motif_list = []

with open(motif_input, "r") as mo:
    for line in mo:
        motif_list.append(line.replace('\n',''))

########################################################################
########## USE MOTIF CONV. FN() TO MAKE DICTIONARY OF REGEX

motif_regex_dict = {}

for m in range(len(motif_list)):
    motif_regex_dict[motif_list[m]] = motif_conversion(motif_list[m])

# motif_regex_dict has keys with original motifs (with ambiguous IUPAC symbols)
# and values with their corresponding reg. expressions 
# built with motif_conversion() fn previously defined

########################################################################
# Creating read and header lists to work from

current_read_list = []
all_reads = []
all_read_headers = []

with open(fasta_input, "r") as fa:
    for line in fa:
        if ">" not in line:
            current_read_list.append(line.replace('\n',''))
        else:
            all_reads.append(''.join(current_read_list))
            current_read_list = []
            all_read_headers.append(line.replace('\n',''))

    all_reads.append(''.join(current_read_list))

# removing first (empty) read
all_reads.pop(0)

# Get list of just the gene symbols of the headers from the FASTA header lines
header_list = [x.split()[0][1:] for x in all_read_headers]


##########################################################
##########################################################
# OUTPUT TO GRAPHIC AND LOOPING THROUGH ALL READS 
# FOR PROCESSING STARTS BELOW HERE
##########################################################




# Setting up cairo surface
num_reads = len(all_reads)
width, height = 1000, (500 * num_reads) + 200

# Setting up output to write to [file prefix of fasta input].svg
surface = cairo.SVGSurface( ("%s.svg" % (fasta_input)), width, height)
ctx = cairo.Context(surface)

# Adding text to annotate motif output

# Setting color, size and style of font
ctx.set_source_rgb(0, 0, 0)
ctx.set_font_size(15)
ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

# Getting dimensions of font output for centering
(x, y, width, height, dx, dy) = ctx.text_extents("Longer vertical bars indicate start of motifs for discerning ambiguous overlapping sites")

# Positioning and outputting text
ctx.move_to( (500 - (width / 2)), ( (500 * num_reads) + 100 - height))
ctx.show_text("Longer vertical bars indicate start of motifs for discerning ambiguous overlapping sites")


# Looping through all the FASTA records, using variable rec (for record)
for rec in range(len(all_reads)):
    aread = all_reads[rec]
    aheader = all_read_headers[rec]

    


    ###############################################################
    # Getting indices from reads - where motifs start
    # Adding to empty (temporary) dictionary

    aread_motif_indices = {}
    
    # Looping through motif_list to check through each motif - adding all to temporary dictionary
    for x in range(len(motif_list)):
        current_motif_indices = [match.start() for match in re.finditer(motif_regex_dict[motif_list[x]], aread)]
        aread_motif_indices[motif_list[x]] = current_motif_indices

    # Above outputs "aread_motif_indices" dictionary with lists as indices of original motifs
    # empty lists indicate that the motif was not found in the read at all

    ###################################################################
    ########### CREATING GRAPHIC WITH PYCAIRO BELOW
    ########### Starting with defining functions
    #
    # Note:
    ## All positions and sizes in graphic output are scaled 
    # to the size of the graphic output, in proportion to the actual genetic information




    x = 0
    #####################################################3
    # Defining color functions
    def intron_exon_col():
        ctx.set_source_rgb(0.541176470588235, 0.772549019607843, 1.0)

    def red():
        ctx.set_source_rgb((247/255), (35/255), (35/255))

    def orange():
        ctx.set_source_rgb((255/255), (140/255), 0)

    def yellow():
        ctx.set_source_rgb((255/255), (234/255), 0)

    def green():
        ctx.set_source_rgb((77/255), (255/255), 0)

    def indigo():
        ctx.set_source_rgb((85/255), (43/255), (117/255))

    def violet():
        ctx.set_source_rgb((252/255), (50/255), (209/255))

    def dark_green():
        ctx.set_source_rgb((46/255), (102/255), (62/255))

    def navy_blue():
        ctx.set_source_rgb((10/255), (92/255), (170/255))

    def burgundy():
        ctx.set_source_rgb((128/255), (60/255), (60/255))

    def purple():
        ctx.set_source_rgb((158/255), (21/255), (142/255))

    def black():
        ctx.set_source_rgb(0, 0, 0)
    ######################################################
    ### Setting up combination function to iterate through color-picking
    def col_fn_list(x):
        if x == 0:
            intron_exon_col()
        
        elif x == 1:
            red()
        
        elif x == 2:
            orange()
        
        elif x == 3:
            yellow()
        
        elif x == 4:
            green()
        
        elif x == 5:
            indigo()
        
        elif x == 6:
            violet()
        
        elif x == 7:
            dark_green()
        
        elif x == 8:
            navy_blue()
        
        elif x == 9:
            burgundy()
        
        elif x == 10:
            purple()
        
        elif x == 11:
            black()
    #########################################################
    # col_fn_list() function chooses color based on number - pulls from previously defined color functions
    # 0: for intron/exon
    # 1 - 10 for colors of motifs
    # 11: black

    ######################################################
    ### Drawing intron to overlay exon on

    # Variables for placement
    # Intron dimensions
    intron_tl_x = 0
    intron_tl_y = 150 + (rec * 500)
    intron_width = 1000
    intron_height = 10


    # Drawing intron
    col_fn_list(0)
    ctx.rectangle(intron_tl_x, intron_tl_y, intron_width, intron_height)
    ctx.fill()

    # Getting positions for exon
    # test_exon_placer_output is a list containing the following: exon sequence, start and end (order-specific)
    test_exon_placer_output = (exon_placer(aread))
    start_pct = (test_exon_placer_output[1] / len(aread)) * 1000
    end_pct = ( (test_exon_placer_output[2] - test_exon_placer_output[1]) / len(aread) ) * 1000


    # Exon
    # Variables for placement
    y_exon = 135 + (rec * 500)
    vert_len_exon = 40

    # Drawing exon
    col_fn_list(0)
    ctx.rectangle(start_pct, y_exon, end_pct, vert_len_exon)
    ctx.fill()
    ######################################################
    ### Adding header of read to top of output

    # Variables for placement
    # header_x_placement is picked based of text_extents, to center text in output
    header_y_placement = 25 + (rec * 500)

    # Making header text - only need gene name - 
    # splitting on whitespace and leaving out the leading ">" character
    gene_name = aheader.split()[0][1:]

    # Setting color, size, style of font
    col_fn_list(11)
    ctx.set_font_size(20)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    # Getting dimensions of text in order to place with respect to center in .svg graphic output
    (x, y, width, height, dx, dy) = ctx.text_extents("Gene Symbol: " + gene_name)

    # Placing and outputting text
    ctx.move_to( (500 - (width / 2)), header_y_placement)
    ctx.show_text("Gene Symbol: " + gene_name)


    ######################################################
    ### Making legend of different motif colors (max of 10)

    ### Building working_motif_list
    # Array containing:
    # list of motif, and corresponding list (in 2nd dimension) of position of occurrences
    # built from previous dictionary - 
    # need list instead of dict to keep ordered (varying versions of python are different in this regard)

    working_motif_list = []

    for (k, v) in aread_motif_indices.items():
        if len(v) > 0:
            alist = [k, v]
            working_motif_list.append(alist)
            


    #####################################################################
    # Making motif legend

    # Variables for placement
    legend_text_x = 50
    legend_text_y = 280 + (rec * 500)

    # Color, size, style of font
    col_fn_list(11)
    ctx.set_font_size(15)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    # Positioning and outputting font
    ctx.move_to(legend_text_x, legend_text_y)
    ctx.show_text("Motif Legend")

    # Making motif legend outline
    # based on length of working_motif_list 
    # how many motifs are found in this sequence

    # Setting dimensions of outline
    left = 50
    right = 175
    top = 290 + (rec * 500)
    bottom = top + (len(working_motif_list) * 25 )


    # Setting color, line width, choosing path and drawing to .svg output
    col_fn_list(11)
    ctx.set_line_width(1)
    ctx.move_to(left, top)
    ctx.line_to(right, top)
    ctx.line_to(right, bottom)
    ctx.line_to(left, bottom)
    ctx.line_to(left, top)
    ctx.close_path()
    ctx.stroke()

    ###########################

    # Setting variable x for further use positioning individual genes in output graphic
    x = 0 + (rec * 500)

    # Setting loop - variable "m" for: 1 motif, and associated information
    for m in range(len(working_motif_list)):
        # Setting color; moving through list of colors in pre-defined swatch
        col_fn_list(m + 1)
        motif_len = len(working_motif_list[m][0])
        curr_motif_width = (motif_len / len(aread)) * 1000

        # Inner loop to move through each of the occurrences of finding a motif
        for p in range(len(working_motif_list[m][1])):
            
            # Setting position for motif start and scaling to size of graphic output
            raw_posn = working_motif_list[m][1][p]
            curr_motif_posn = (raw_posn / len(aread) ) * 1000

            # Setting dimensions for motif
            y_top_main = 130 + (rec * 500)
            y_main_length = 50
            y_top_starting = 100 + (rec * 500)
            y_starting_length = 110


            # Drawing motifs, composed of two parts:
            # Main portion - the entire length of the motif is represented by this rectangle
            ctx.rectangle(curr_motif_posn, y_top_main, curr_motif_width, y_main_length)
            # "Starting" portion - the exact start of the motif is marked by this thin "line" (also a rectangle)
            ctx.rectangle(curr_motif_posn, y_top_starting, 1, y_starting_length)

        # Making swatches of color (squares) in the key legend, 
        # to annotate motif colors with sequences from input

        # Variables for placement
        legend_col_x = 60

        # Drawing squares to output
        ctx.rectangle(legend_col_x, (300 + x), 10, 10)
        ctx.fill()

        # Annotating legend with motif sequences:

        # Variables for placement
        motif_seqs_x = 85

        # Choosing color, size and style of font
        col_fn_list(11)
        ctx.set_font_size(10)
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

        # Writing out individual motifs to their corresponding color swatches in the legend
        # Moving through system in a loop
        ctx.move_to(motif_seqs_x, (310 + x) )
        ctx.show_text(working_motif_list[m][0])


        # Incrementing x (used in placement)
        x += 15